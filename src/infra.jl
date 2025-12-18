"""
AMM base interface. An AMM has parameters applied but has not been initialized at a given liquidity.
Corresponds to the "CFMM" concept in the liquidity density paper.

Expected Interface:
- `kind_name(self)::str`
- `short_str(self, as_bp::bool)::str`
TODO ^^^
"""
abstract type AMM end

"""
    kind_name(amm)
    kind_name(tc)

Name of the AMM type (without parameters)

Optional.
"""
kind_name(::T) where {T<:AMM} = string(nameof(T))

"""
    display_params(::AMM; as_bp::Bool)

Must return a list of pairs of strings, for display only. Used in the default implementation of `short_str()`.

Required.
"""
display_params(::AMM; as_bp::Bool) = []

# Helper
render_display_params(params) = join(["$k=$v" for (k, v) in params], ", ")

"""
    short_str(amm; as_bp::Bool)
    short_str(tc; as_bp::Bool)

Short description of the AMM with parameters. If `as_bp=true`, render all prices in bp.

Optional, rarely needs an overload
"""
function short_str(amm::AMM; as_bp::Bool)
    # SOMEDAY I think this should actually asked for another method that returns a list of params in a simpler way. That other method can then be overwritten.
    params = display_params(amm; as_bp = as_bp)
    params_str = render_display_params(params)
    name = kind_name(amm)
    "$name($params_str)"
end

"""
    alpha, beta = alphabeta(amm)
    alpha, beta = alphabeta(tc)

The lower and upper end of the price range, or (0.0, Inf).

Required for AMMs with concentrated price range.
"""
alphabeta(::AMM) = (0.0, Inf)
# The defaults are full range. AMMs need to override this if not.

"""
    t_l_p(::AMM, l, p)

Balances at the given invariant l and price p. 

Required.
"""
function t_l_p(::AMM, l, p)
    error("Abstract base method")
end

"""
    p_l_t(::AMM, l, t)

Price at the given invariant l and balances t.

Required.
"""
function p_l_t(::AMM, l, t)
    error("Abstract base method")
end

"""
    l_t(::AMM, t)

Invariant at the given balances t.

Required.
"""
function l_t(::AMM, t) end

# TODO implement this for other functions. (I think this is done??)
# Should be differentiable if possible.
# NB The default implementation could be `l_t(t) - l` but this is usually not efficient.
# This is _not_ non-negative. Square it or use abs if you need that.
"""
    residual(amm, l, t)
    residual(tc, t)

Invariant residual at the given invariant l and balances t. This is 0 iff t is on the given trading curve.

TODO Required? Should be differentiable.
"""
function residual end

"""
    tc_t(amm::AMM, t) <: TradingCurve

Trading curve of an AMM at the given balances t.
"""
function tc_t(amm::AMM, t)
    l = l_t(amm, t)
    p0 = p_l_t(amm, l, t)
    mk_tc(TCCommon(amm, l, t, p0))
end

"""
    tc_v_p(amm::AMM, v, p)

Trading curve of an AMM at the given portfolio value v and price p.
"""
function tc_v_p(amm::AMM, v, p0)
    l, t = lt_vp(amm, v, p0)
    mk_tc(TCCommon(amm, l, t, p0))
end

# NB the trading curve itself doesn't need the price but it has a p0 member.
# Maybe it shouldn't or it should be optional, idk.
"""
    tc_l_p(amm::AMM, l, p)

Trading curve of an AMM at the given invariant l and price p
"""
function tc_l_p(amm::AMM, l, p0)
    t = t_l_p(amm, l, p0)
    mk_tc(TCCommon(amm, l, t, p0))
end

# backwards-compatible aliases
const trading_curve_from_balances = tc_t
const trading_curve_from_value_price = tc_v_p
const trading_curve_from_invariant_price = tc_l_p

"""
  flip(amm)
  flip(tc)

Flip the x and y assets.

Flipping an AMM changes what the two assets mean. Flipping a TradingCurve yields the _same_ TradingCurve with the same initial state but with the x and y assets flipped.
"""
function flip end

"""
    lt_vp(amm::AMM, v, p)

Compute the invariant and balances at a given portfolio value (in units of the y asset) and price.

The fallback implementation requires `t_l_p` and scales the invariant from L=1.
"""
function lt_vp(amm::AMM, v, p)
    # TODO check against price bounds

    t1 = t_l_p(amm, 1.0, p)
    pfvaly1 = dot((p, 1.0), t1)

    l = v / pfvaly1
    return l, t1 .* l
end

export TradingCurve, TCCommon, t_p, x_y, y_x, t_p, dydp_p

"""
A trading curve. It is expected that this has a field `.common` of type `TCCommon`.

SOMEDAY that's kinda bad actually.
"""
abstract type TradingCurve end

# Make AMM and TradingCurve behave like a scalar in broadcasts (useful for syntax)
Base.Broadcast.broadcastable(amm::AMM) = Ref(amm)
Base.Broadcast.broadcastable(tc::TradingCurve) = Ref(tc)

"""
    t_plus(::TradingCurve)

Maximal x and and y values along the trading curve. Inf if none, respectively.

TODO refactor: this gives an error if not implemented but alphabeta defaults to (0, Inf). Inconsistent. This behavior is prob better.
"""
t_plus(::TradingCurve) = error("Abstract base method")

kind_name(tc::TradingCurve) = kind_name(tc.common.amm)

"""
    t_p(tc::TradingCurve, p)

Balances at given price p
"""
t_p(tc::TradingCurve, p) = t_l_p(tc.common.amm, tc.common.l, p)

alphabeta(tc::TradingCurve) = alphabeta(tc.common.amm)

"""
    x_y(tc, y)

x balance given the y balance.
"""
function x_y end

"""
    y_x

y balance given the x balance.
"""
function y_x end

struct TCCommon{A<:AMM}
    amm::A
    l::Float64
    t_init::Vector{Float64}
    p_init::Float64
end

"""
    v_init(common::TCCommon)
    v_init(tc)

Portfolio value / TVL denominated in y asset units, at pool initialization
"""
v_init(common::TCCommon) = dot(common.t_init, [common.p_init, 1.0])
v_init(tc::TradingCurve) = v_init(tc.common)

function short_str(tc::TradingCurve; as_bp::Bool)
    amm_str = short_str(tc.common.amm; as_bp)
    pfval_str = fmt_auto(v_init(tc))
    "$amm_str @ \$$pfval_str"
end

"""
    mk_tc(common) <: TradingCurve

Create a certain trading curve object from an AMM and common data. Must be specialized to the concrete type of `common`. This is how we capture the relationship from AMMs to TradingCurve's

Required.
"""
mk_tc(::TCCommon) = error("Abstract base method.")

"""
    dydp_p(tc, p)

Absolute liquidity density at the given price

NOTE: This is only defined for AMMs where this is fast to compute. Don't use this for plotting; instead, use `sample_dydp_p()`.
"""
function dydp_p end

"""
    dydlogp_p(tc, p)

Relative liquidity density at the given price

NOTE: This is only defined for AMMs where this is fast to compute. Don't use this for plotting; instead, use `sample_dydlogp_p()`.
"""
function dydlogp_p end

"""
    sample_dydp_p(tc)

Sample (p, dydp) points (absolute LD curve) from the given TradingCurve. Returns a vector of pairs.

# Keyword parameters

- `plim`: If given, limit to that range. *must* be given for trading curves with infinite price range.
- `unit`: Unit of price change. The default is 1.0 = ±1.
"""
function sample_dydp_p(tc::TradingCurve; plim = (0.0, Inf), unit = 1.0)
    # TODO mogrify this into returning nothing for empty, like the other instance of the function, and move that function actually.
    plim = Util.intersect_intervals(alphabeta(tc), plim)
    @assert plim[1] < plim[2]
    @assert plim[1] > 0.0
    @assert isfinite(plim[2])
    res = _sample_dydp_p(tc, plim[1], plim[2])
    map(res) do (p, v)
        (p, v * unit)
    end
end

"""
    sample_dydlogp_p(tc)

Sample (p, dydlogp) points (relative LD curve) from the given TradingCurve. Returns a vector of pairs.

# Keyword parameters

- `plim`: If given, limit to that range. *must* be given for trading curves with infinite price range.
- `unit`: Unit of price change. The default is 0.01 = ±1%.
"""
function sample_dydlogp_p(tc::TradingCurve; plim = (0.0, Inf), unit = 0.01)
    plim = Util.intersect_intervals(alphabeta(tc), plim)
    @assert plim[1] < plim[2]
    @assert plim[1] > 0.0
    @assert isfinite(plim[2])
    res = _sample_dydlogp_p(tc, plim[1], plim[2])
    map(res) do (p, v)
        (p, v * log(1.0 + unit))
    end
end

# Unchecked, untransformed versions for overloading
# The default implementation uses dydp_p. Must be overwritten for others.
_sample_dydp_p(tc::TradingCurve, alpha, beta) =
    sample_adaptive(p -> dydp_p(tc, p), alpha, beta)
_sample_dydlogp_p(tc::TradingCurve, alpha, beta) =
    sample_adaptive(p -> dydp_p(tc, p) * p, alpha, beta)

"""
    dydp_max(tc; plim=(0.0, Inf)) :: Tuple{Float64,Float64}

(price, value) of the maximal value of dydp (absolute liquidity density), optionally within the
interval `plim`. If `plim` is not a subset of the price set, the intersection is taken.

# Keyword parameters

- `plim`: If given, limit to that range. *must* be given for trading curves with infinite price range.
- `unit`: Unit of price change. The default is 1.0 = ±1.
- `pathological`: If true, we include pathological limit behavior for some AMMs (e.g., StableSwap's
  infinite LD if plim includes 0.0). Otherwise, these points are ignored.
"""
function dydp_max(
    tc::TradingCurve;
    plim = (0.0, Inf),
    strict = false,
    unit = 1.0,
    pathological = false,
)
    plim = Util.intersect_intervals(alphabeta(tc), plim)
    @assert plim[1] < plim[2]
    (p, v) = _dydp_max(tc, plim...; pathological)
    v *= unit
    (p, v)
end

"""
    dydlogp_max(tc::TradingCurve; plim=(0.0, Inf)) :: Tuple{Float64,Float64}

(price, value) of the value of dydlogp (relative liquidity density), optionally within the interval
`plim`. If `plim` is not a subset of the price set, the intersection is taken.

# Keyword parameters

- `plim`: If given, limit to that range. *must* be given for trading curves with infinite price range.
- `unit`: Unit of price change. The default is ±1%.
- `pathological`: If true, we include pathological limit behavior for some AMMs (e.g., StableSwap's
  infinite LD if plim includes 0.0). Otherwise, these points are ignored.
"""
function dydlogp_max(tc::TradingCurve; plim = (0.0, Inf), unit = 0.01, pathological = false)
    plim = Util.intersect_intervals(alphabeta(tc), plim)
    @assert plim[1] < plim[2]
    (p, v) = _dydlogp_max(tc, plim...; pathological)
    v *= log(1.0 + unit)
    (p, v)
end

using Optim

# The default implementations simply use numeric optimization using dydp_p and ignores the
# `pathological`. This may be suboptimal or illegal in some instances, then this needs to be
# overwritten, in particular for infinite price ranges.
function _dydp_max(tc::TradingCurve, alpha, beta; pathological)
    # NB optimize always computes the minimum.
    res = optimize(p -> -dydp_p(tc, p), alpha, beta)
    @assert res.stopped_by.converged
    (res.minimizer, -res.minimum)
end

function _dydlogp_max(tc::TradingCurve, alpha, beta; pathological)
    # NB optimize always computes the minimum.
    res = optimize(p -> -dydp_p(tc, p) * p, alpha, beta)
    @assert res.stopped_by.converged
    (res.minimizer, -res.minimum)
end

flip(tc::TradingCurve) = tc_t(flip(tc.common.amm), reverse(tc.common.t_init))

residual(tc::TradingCurve, t) = residual(tc.common.amm, tc.common.l, t)
