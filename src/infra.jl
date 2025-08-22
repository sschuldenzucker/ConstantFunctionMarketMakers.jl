"""
AMM base interface. An AMM has parameters applied but has not been initialized at a given liquidity.
Corresponds to the "CFMM" concept in the liquidity density paper.

Expected Interface:
- `kind_name(self)::str`
- `short_str(self, as_bp::bool)::str`
"""
abstract type AMM end

"""
Name of the AMM type (without parameters)
"""
kind_name(::T) where {T<:AMM} = string(nameof(T))

"""
Must return a list of pairs of strings, for display only. Used in the default implementation of `short_str()`.
"""
display_params(::AMM; as_bp::Bool) = []

# Helper
render_display_params(params) = join(["$k=$v" for (k, v) in params], ", ")

"""
Short description of the AMM with parameters. If `as_bp=True`, render all prices in bp.
"""
# SOMEDAY I think this should actually asked for another method that returns a list of params in a simpler way. That other method can then be overwritten.
function short_str(amm::AMM; as_bp::Bool)
    params = display_params(amm; as_bp = as_bp)
    params_str = render_display_params(params)
    name = kind_name(amm)
    "$name($params_str)"
end

"""
The lower and upper end of the price range, or (0.0, Inf)
"""
function alphabeta end

# The defaults are full range. AMMs need to override this if not.
alphabeta(::AMM) = (0.0, Inf)

"""
`t_l_p(<:AMM, l, p)` returns the balances at the given invariant l and price p. 
"""
function t_l_p end

function p_l_t end

function l_t end

# TODO implement this for other functions.
# Should be differentiable if possible.
# NB The default implementation could be `l_t(t) - l` but this is usually not efficient.
# This is _not_ non-negative. Square it or use abs if you need that.
function residual end

function tc_t(amm::AMM, t)
    l = l_t(amm, t)
    p0 = p_l_t(amm, l, t)
    mk_tc(TCCommon(amm, l, t, p0))
end

function tc_v_p(amm::AMM, v, p0)
    l, t = lt_vp(amm, v, p0)
    mk_tc(TCCommon(amm, l, t, p0))
end

# NB the trading curve itself doesn't need the price but it has a p0 member.
# Maybe it shouldn't or it should be optional, idk.
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

function t_plus end

kind_name(tc::TradingCurve) = kind_name(tc.common.amm)

t_p(tc::TradingCurve, p) = t_l_p(tc.common.amm, tc.common.l, p)

alphabeta(tc::TradingCurve) = alphabeta(tc.common.amm)

function x_y end
function y_x end

struct TCCommon{A<:AMM}
    amm::A
    l::Float64
    t_init::Vector{Float64}
    p_init::Float64
end

v_init(common::TCCommon) = dot(common.t_init, [common.p_init, 1.0])
v_init(tc::TradingCurve) = v_init(tc.common)

function short_str(tc::TradingCurve; as_bp::Bool)
    amm_str = short_str(tc.common.amm; as_bp)
    pfval_str = fmt_auto(v_init(tc))
    "$amm_str @ \$$pfval_str"
end

"""
    mk_tc(common, ...) <: TradingCurve

Trait-like function to create a certain trading curve object from an AMM and common data. Must be specialized to the concrete type of `common` 

SOMEDAY this doesn't actually need its `amm` parameter anymore b/c common.amm stores it and it's type-specialized now.
"""
mk_tc(::TCCommon) = error("Abstract base method.")

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
    dydp_max(tc::TradingCurve; plim=(0.0, Inf)) :: Tuple{Float64,Float64}

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
