
export concentrate, ConcentratedAMM, ConcentratedTradingCurve

# TODO inconsistent how this works for trading curves but e.g. rate_scale does not. Or is this the only one where this makes actual sense?
"""
    concentrate(amm, alpha, beta) <: AMM
    concentrate(tc, alpha, beta)
    amm |> concentrate(alpha, beta)
    tc |> concentrate(alpha, beta)

Concentrate the given AMM or TradingCurve into a sub-range. This corresponds to (TOOD Theorem ref). `(alpha, beta)` must be a sub-interval of the price range.

For AMMs, returns a new AMM concentrated into that sub-range.

For trading curves, returns the corresponding curve when preserving the given trading curve's initial price and value. (the initial price must be within the new price bounds)
"""
function concentrate(amm::AMM, alpha, beta)
    @assert Util.is_subinterval((alpha, beta), alphabeta(amm))
    _concentrate(amm, alpha, beta)
end

function concentrate(tc::TradingCurve, alpha, beta)
    @assert Util.is_subinterval((alpha, beta), alphabeta(tc))
    p_init = tc.common.p_init
    @assert alpha <= p_init <= beta "TC init price must be within new price bounds"
    tc_v_p(concentrate(tc.common.amm, alpha, beta), v_init(tc), p_init)
end

# piping helper
concentrate(alpha, beta) = x -> concentrate(x, alpha, beta)

_concentrate(::CPMM, alpha1, beta1) = C2LP(alpha1, beta1)

# TODO check if these are true. I think they are Use my theorem for concentrating any AMM.
# There seems to be an invariant here. ("concentrating twice" or something)

_concentrate(::C2LP, alpha1, beta1) = C2LP(alpha1, beta1)  # TODO check this but should be true

function _concentrate(amm::ECLP, alpha1, beta1)
    # This is not strictly required but I _think_ we may need it for our implementation, and we also check it in the ECLP constructor. I think
    @assert alpha1 <= amm.params.q <= beta1
    ECLP(alpha1, beta1, amm.params.q, amm.params.lam)
end

# The default implementation uses ConcentratedAMM.
_concentrate(amm::AMM, alpha, beta) = ConcentratedAMM(amm, alpha, beta)

"""
Generic concentrated AMM, using numerical methods to solve the required equations.

Note: `concentrate()` uses knowledge about some AMMs to yield more efficient versions.
`ConcentratedAMM()` always uses the AMM as an approximate black box. (though some
optimization/parameter choice may depend on the inner AMM)
"""
struct ConcentratedAMM{A<:AMM} <: AMM
    inner::A
    alpha::Float64
    beta::Float64

    chi::Vector{Float64}

    function ConcentratedAMM(inner::A, alpha::Float64, beta::Float64) where {A<:AMM}
        # NB chi is always finite by construction even if alpha or beta are infinite.
        chi = [t_l_p(inner, 1.0, beta)[1], t_l_p(inner, 1.0, alpha)[2]]
        new{A}(inner, alpha, beta, chi)
    end
end

flip(amm::ConcentratedAMM) = ConcentratedAMM(flip(amm.inner), 1 / amm.beta, 1 / amm.alpha)

kind_name(amm::ConcentratedAMM) = "Concentrated($(kind_name(amm.inner)))"

function short_str(amm::ConcentratedAMM; as_bp)
    inner_str = short_str(amm.inner; as_bp = as_bp)
    alpha_s, beta_s = fmt_auto.((amm.alpha, amm.beta))
    "$inner_str ∩ ($alpha_s, $beta_s)"
end

alphabeta(amm::ConcentratedAMM) = (amm.alpha, amm.beta)

residual(amm::ConcentratedAMM, l, t) = residual(amm.inner, @. t + l * amm.chi)

function l_t(amm::ConcentratedAMM, t)
    # NB this default implementation just works for all AMMs we have right now but may need adjustments to the method for specific AMMs. In that case, we'll want to specialize further.
    # SOMEDAY the starting point may be bad. Might instead use `l_t(amm.inner, t)` or something custom.
    # Also note that we could solve most of these analytically.
    res = nlsolve(v -> residual(amm, v[1], t), [1.0]; autodiff = :forward)
    @assert res.f_converged
    res.zero[1]
end

t_l_p(amm::ConcentratedAMM, l, p) = t_l_p(amm.inner, l, p) .- l .* amm.chi
p_l_t(amm::ConcentratedAMM, l, t) = p_l_t(amm.inner, l, @. t + l * amm.chi)

function mk_tc(common::TCCommon{ConcentratedAMM{A}}) where {A}
    amm = common.amm
    common_inner = TCCommon(
        amm.inner,
        common.l,
        (@. common.t_init + common.l * amm.chi),
        common.p_init,
    )
    ConcentratedTradingCurve(common, mk_tc(common_inner))
end

struct ConcentratedTradingCurve{A<:AMM,T<:TradingCurve} <: TradingCurve
    common::TCCommon{ConcentratedAMM{A}}
    inner::T
end

t_plus(tc::ConcentratedTradingCurve) =
    let chi = tc.common.amm.chi, l = tc.common.l
        t_plus(tc.inner) .- l .* chi
    end

x_y(tc::ConcentratedTradingCurve, y) =
    let chi = tc.common.amm.chi, l = tc.common.l
        x_y(tc.inner, y + l * chi[2]) - l * chi[1]
    end

y_x(tc::ConcentratedTradingCurve, x) =
    let chi = tc.common.amm.chi, l = tc.common.l
        x_y(tc.inner, x + l * chi[1]) - l * chi[2]
    end

# NB only defined if dydp_p is defined for inner.
dydp_p(tc::ConcentratedTradingCurve, p) = dydp_p(tc.inner, p)
dydp_y(tc::ConcentratedTradingCurve, y) = dydp_y(tc.inner, y)

_sample_dydp_p(tc::ConcentratedTradingCurve, alpha, beta) =
    _sample_dydp_p(tc.inner, alpha, beta)
_sample_dydlogp_p(tc::ConcentratedTradingCurve, alpha, beta) =
    _sample_dydlogp_p(tc.inner, alpha, beta)
_dydp_max(tc::ConcentratedTradingCurve, alpha, beta; pathological) =
    _dydp_max(tc.inner, alpha, beta; pathological)
_dydlogp_max(tc::ConcentratedTradingCurve, alpha, beta; pathological) =
    _dydlogp_max(tc.inner, alpha, beta; pathological)

