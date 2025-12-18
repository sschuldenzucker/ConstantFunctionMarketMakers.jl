export rate_scale, RateScaledAMM, RateScaledTradingCurve, constant_rate_scaled_eclp

# Convenience helper for piping
# SOMEDAY make it operate on TCs, too? Cf. concentrate(). Not super clear what the right basis is, maybe t_init??
"""
    rate_scale(amm, rx, ry)
    amm |> rate_scale(rx, ry)

Return a new AMM, which is the given AMM scaled by the given rates.
"""
rate_scale(amm::AMM, rx, ry) = RateScaledAMM(amm, rx, ry)
rate_scale(rx, ry) = x -> rate_scale(x, rx, ry)

# SOMEDAY also the number types should maybe be flexible.
# To calculate derivates, for example.
struct RateScaledAMM{A<:AMM} <: AMM
    inner::A
    # Someday should these be scaling_factors::Vector{Float64}?
    rx::Float64
    ry::Float64
end

flip(amm::RateScaledAMM) = RateScaledAMM(flip(amm.inner), amm.ry, amm.rx)

# NB we don't use display_params or kind_name for display.

function kind_name(rsa::RateScaledAMM)
    inner_name = kind_name(rsa.inner)
    "Scaled($inner_name)"
end

function short_str(rsa::RateScaledAMM; as_bp)
    inner_str = short_str(rsa.inner; as_bp = as_bp)
    srx, sry = fmt_auto.((rsa.rx, rsa.ry))
    "$inner_str * ($srx, $sry)"
end

# Specialization to get nice-looking ECLP labels
# NB this is a bit of a hack.
function short_str(rsa::RateScaledAMM{ECLP}; as_bp)
    # Still wondering if display_params() should just get a rate scaling parameter instead or so.
    eclp = rsa.inner
    mogrify = as_bp ? fmt_todev : fmt_auto
    params = [
        "α" => scale_inner_price(rsa, eclp.params.alpha) |> mogrify,
        "β" => scale_inner_price(rsa, eclp.params.beta) |> mogrify,
        "peg" => scale_inner_price(rsa, eclp.params.q) |> mogrify,
        "λ" => eclp.params.lam,
    ]
    params_str = render_display_params(params)
    name = kind_name(eclp)
    "$name($params_str) (scaled)"
end

scale_inner_price(rsa::RateScaledAMM, p) = rsa.rx / rsa.ry * p
scale_outer_price(rsa::RateScaledAMM, p) = rsa.ry / rsa.rx * p
scale_inner_balances(rsa::RateScaledAMM, t) = t ./ [rsa.rx, rsa.ry]
scale_outer_balances(rsa::RateScaledAMM, t) = [rsa.rx, rsa.ry] .* t

alphabeta(rsa::RateScaledAMM) = scale_inner_price.(rsa, alphabeta(rsa.inner))

# SOMEDAY technically the inner functions are not always defined (e.g. for StableSwap). We could use traits to check for this and throw a method error earlier but I can't be bothered right now.

residual(rsa::RateScaledAMM, l, t) = residual(rsa.inner, l, scale_outer_balances(rsa, t))

l_t(rsa::RateScaledAMM, t) = l_t(rsa.inner, scale_outer_balances(rsa, t))

t_l_p(rsa::RateScaledAMM, l, p) =
    scale_inner_balances(rsa, t_l_p(rsa.inner, l, scale_outer_price(rsa, p)))

p_l_t(rsa::RateScaledAMM, l, t) =
    scale_inner_price(rsa, p_l_t(rsa.inner, l, scale_outer_balances(rsa, t)))

function mk_tc(common::TCCommon{RateScaledAMM{A}}) where {A}
    # NB this is correct but a bit 'von hinten durch die Brust ins Auge' b/c we have previously
    # scaled the inner balances to the outer balances order, and now we're scaling back. We could
    # just specialize the two `trading_curve_from_*` functions instead.
    rsa = common.amm
    common_inner = TCCommon(
        rsa.inner,
        common.l,
        scale_outer_balances(rsa, common.t_init),
        scale_outer_price(rsa, common.p_init),
    )
    RateScaledTradingCurve(common, mk_tc(common_inner))
end

struct RateScaledTradingCurve{A<:AMM,T<:TradingCurve} <: TradingCurve
    common::TCCommon{RateScaledAMM{A}}
    # Note that the inner AMM is duplicated, once as inner.common.amm, and once as common.amm.inner.
    inner::T
end

t_plus(tc::RateScaledTradingCurve) = scale_inner_balances(tc.common.amm, t_plus(tc.inner))

x_y(tc::RateScaledTradingCurve, y) = x_y(tc.inner, y * tc.common.amm.ry) / tc.common.amm.rx
y_x(tc::RateScaledTradingCurve, x) = y_x(tc.inner, x * tc.common.amm.rx) / tc.common.amm.ry

# NOTE this will only succeed if inner has t_p defined (i.e., not for StableSwap)
# The same is true for some of the methods above.
# TODO is there some better way. Trait checks don't work here either. Maybe manual trait?
t_p(tc::RateScaledTradingCurve, p) =
    scale_inner_balances(tc.common.amm, t_p(tc.inner, scale_outer_price(tc.common.amm, p)))

# TODO Same issue as for t_p.
# Yes, this formula is asymmetric!!
dydp_p(tc::RateScaledTradingCurve, p) =
    dydp_p(tc.inner, scale_outer_price(tc.common.amm, p)) / tc.common.amm.rx

# Only if dydp_y is defined for the underlying AMM.
# NB this can have a default implementation but is only ever needed as a fallback for StableSwap.
dydp_y(tc::RateScaledTradingCurve, y) =
    dydp_y(tc.inner, y * tc.common.amm.ry) / tc.common.amm.rx

function _sample_dydp_p(tc::RateScaledTradingCurve, alpha, beta)
    # SOMEDAY in principle, sampling the inner one should have different parameters (max error etc.) than the outer one.
    amm = tc.common.amm
    points_inner = _sample_dydp_p(
        tc.inner,
        scale_outer_price(amm, alpha),
        scale_outer_price(amm, beta),
    )
    [(scale_inner_price(amm, p), dydp / amm.rx) for (p, dydp) in points_inner]
end

function _sample_dydlogp_p(tc::RateScaledTradingCurve, alpha, beta)
    # SOMEDAY in principle, sampling the inner one should have different parameters (max error etc.) than the outer one.
    amm = tc.common.amm
    points_inner = _sample_dydlogp_p(
        tc.inner,
        scale_outer_price(amm, alpha),
        scale_outer_price(amm, beta),
    )
    [
        # Note that a different scaling factor is used here! ry instead of rx!
        (scale_inner_price(amm, p), dydp / amm.ry) for (p, dydp) in points_inner
    ]
end

function _dydp_max(tc::RateScaledTradingCurve, alpha, beta; pathological)
    amm = tc.common.amm
    (p_inner, v) = _dydp_max(
        tc.inner,
        scale_outer_price(amm, alpha),
        scale_outer_price(amm, beta);
        pathological,
    )
    (scale_inner_price(amm, p_inner), v / amm.rx)
end

function _dydlogp_max(tc::RateScaledTradingCurve, alpha, beta; pathological)
    amm = tc.common.amm
    (p_inner, v) = _dydlogp_max(
        tc.inner,
        scale_outer_price(amm, alpha),
        scale_outer_price(amm, beta);
        pathological,
    )
    (scale_inner_price(amm, p_inner), v / amm.ry)
end

# Application: Rate-scaled ECLP
"""
    constant_rate_scaled_eclp(alpha, beta, q; q_inner_override)  

Define a rate-scaled ECLP scaled by a constant rate. `alpha`, `beta`, and `q` refer to "outer"
prices. The inner parameters will be chosen such that these outer prices are exhibited, and
the inner `q` is 1.0. Set `q_inner` to change that. This can be used to adjust the peak (since it's
not exactly at q, only approximately). `lam` always refers to the inner curve.
"""
function constant_rate_scaled_eclp(alpha, beta, q, lam; q_inner = 1.0)
    alpha_inner = alpha / q
    beta_inner = beta / q
    RateScaledAMM(ECLP(alpha_inner, beta_inner, q_inner, lam), q, 1.0)
end

