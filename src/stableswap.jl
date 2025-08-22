include("./stableswap_math.jl")
import .StableSwapMath

export StableSwap, StableSwapTradingCurve

# TODO rename AMM to CFMM to match the paper.
"""
NOTE: Only 2 assets are supported right now and also assumed.
"""
struct StableSwap <: AMM
    a::Float64

    # We convert between the "code style" of A (used externally) and the "math style" of A (used
    # internally). SOMEDAY that is a bit of a mess and could easily be fixed.
    # TODO this might actually be pretty broken.
    StableSwap(a) = new(a / 2)
end

flip(amm::StableSwap) = amm

function display_params(ss::StableSwap; as_bp::Bool)
    ["A" => ss.a |> fmt_auto]
end

alphabeta(::StableSwap) = (0, Inf)

residual(ss::StableSwap, l, t) = StableSwapMath.residual_stableswap(ss.a, t, l)

l_t(ss::StableSwap, t) = StableSwapMath.find_invariant(ss.a, t)

# NB this function is not very fast.
function t_l_p(ss::StableSwap, l, p)
    # cache = StableSwapMath.Cache(ss.a, l)
    # StableSwapMath.find_t_l_p(cache, l, p)
    StableSwapMath.find_t_l_p(ss.a, l, p)
end

function p_l_t(ss::StableSwap, l, t)
    # cache = StableSwapMath.Cache(ss.a, l)
    # cache.p_y(t[2])
    StableSwapMath.p_y_gen(ss.a, l, t[2])
end

function mk_tc(common::TCCommon{StableSwap})
    # cache = StableSwapMath.Cache(common.amm.a, common.l)
    # StableSwapTradingCurve(common, cache)
    StableSwapTradingCurve(common)
end

struct StableSwapTradingCurve <: TradingCurve
    common::TCCommon{StableSwap}
    # _cache::StableSwapMath.Cache
end

x_y(tc::StableSwapTradingCurve, y) = StableSwapMath.x_y_gen(tc.common.amm.a, tc.common.l, y)
y_x(tc::StableSwapTradingCurve, x) = StableSwapMath.y_x_gen(tc.common.amm.a, tc.common.l, x)
# y_x(tc::StableSwapTradingCurve, x) = tc._cache.y_x(x)

# Special for StableSwap
p_y(tc::StableSwapTradingCurve, y) = StableSwapMath.p_y_gen(tc.common.amm.a, tc.common.l, y)
# dydp_y(tc::StableSwapTradingCurve, y) = tc._cache.dydp_y(y)
dydp_y(tc::StableSwapTradingCurve, y) =
    StableSwapMath.dydp_y_gen(tc.common.amm.a, tc.common.l, y)

# NB this function is not very fast. You prob don't wanna use it in a loop but sample y instead.
# SOMEDAY do we need it?
# t_p(tc::StableSwapTradingCurve, p) = StableSwapMath.find_t_l_p(tc._cache, tc.common.l, p)
t_p(tc::StableSwapTradingCurve, p) =
    StableSwapMath.find_t_l_p(tc.common.amm.a, tc.common.l, p)

t_plus(tc::StableSwapTradingCurve) = (Inf, Inf)

function _sample_dydp_p(tc::StableSwapTradingCurve, alpha, beta)
    f(y) = (p_y(tc, y), dydp_y(tc, y))
    y_start = t_p(tc, alpha)[2]
    y_end = t_p(tc, beta)[2]
    _, points = sample_adaptive_parametric(f, y_start, y_end)
    points
end

function _sample_dydlogp_p(tc::StableSwapTradingCurve, alpha, beta)
    f(y) =
        let p = p_y(tc, y)
            (p, dydp_y(tc, y) * p)
        end
    y_start = t_p(tc, alpha)[2]
    y_end = t_p(tc, beta)[2]
    _, points = sample_adaptive_parametric(f, y_start, y_end)
    points
end

function _dydp_max(tc::StableSwapTradingCurve, alpha, beta; pathological)
    # Some work b/c we don't have dydp_p for stableswap and also an infinite price range.
    # SOMEDAY we can exploit the knowledge that stuff is maximized at 1.0 but I'm not sure enough it's exact yet, lol.

    if alpha == 0.0
        if pathological
            return (0.0, Inf)
        end
        alpha = min(1e-4, beta / 2.0)
    end

    if beta == Inf
        beta = max(1e4, alpha * 2.0)
    end

    y_start = t_p(tc, beta)[2]
    y_end = t_p(tc, alpha)[2]

    res = optimize(y -> -dydp_y(tc, y), y_start, y_end)
    @assert res.stopped_by.converged

    (p_y(tc, res.minimizer), -res.minimum)
end

function _dydlogp_max(tc::StableSwapTradingCurve, alpha, beta; pathological)
    if alpha == 0.0
        alpha = min(1e-4, beta / 2.0)
    end

    if beta == Inf
        if pathological
            return (Inf, Inf)
        end
        beta = max(1e4, alpha * 2.0)
    end

    y_start = t_p(tc, beta)[2]
    y_end = t_p(tc, alpha)[2]

    res = optimize(y -> -dydp_y(tc, y), y_start, y_end)
    @assert res.stopped_by.converged

    (p_y(tc, res.minimizer), -res.minimum)
end

