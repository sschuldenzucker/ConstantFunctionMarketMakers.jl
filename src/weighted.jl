
export WeightedPool, WeightedPoolTradingCurve

"""
A weighted product pool (like Balancer's weighted pool), following \$x^delta * y^(1-delta) = L\$.  

This is equivalent to the CPMM at delta = 0.5.

delta must be strictly between 0 and 1.
"""
struct WeightedPool <: AMM
    delta::Float64

    function WeightedPool(delta)
        @assert 0 < delta < 1
        new(delta)
    end
end

flip(amm::WeightedPool) = WeightedPool(1.0 - amm.delta)

display_params(wp::WeightedPool; as_bp::Bool) = ["δ" => fmt_auto(wp.delta)]

alphabeta(::WeightedPool) = (0.0, Inf)

residual(wp::WeightedPool, l, t) = t[1]^wp.delta * t[2]^(1.0 - wp.delta) - l

l_t(wp::WeightedPool, t) = t[1]^wp.delta * t[2]^(1.0 - wp.delta)

function t_l_p(wp::WeightedPool, l, p)
    delta = wp.delta
    deltabar = 1.0 - wp.delta
    theta = (1.0 / (1.0 - wp.delta) - 1.0)
    l .* [p^(-deltabar) * theta^deltabar, p^delta * theta^(-delta)]
end

function p_l_t(wp::WeightedPool, l, t)
    theta = (1.0 / (1.0 - wp.delta) - 1.0)
    theta * t[1] / t[2]
end

mk_tc(common::TCCommon{WeightedPool}) = WeightedPoolTradingCurve(common)

struct WeightedPoolTradingCurve <: TradingCurve
    common::TCCommon{WeightedPool}
end

x_y(tc::WeightedPoolTradingCurve, y) = tc.common.l / y^(1 - tc.common.amm.delta)
y_x(tc::WeightedPoolTradingCurve, x) = tc.common.l / x^tc.common.amm.delta

t_plus(::WeightedPoolTradingCurve) = (Inf, Inf)

dydp_p(tc::WeightedPoolTradingCurve, p) = ForwardDiff.derivative(p -> t_p(tc, p)[2], p)

function _dydp_max(tc::WeightedPoolTradingCurve, alpha, beta; pathological)
    if alpha == 0.0
        return Inf
    end
    (alpha, dydp_p(tc, alpha))
end

function _dydlogp_max(tc::WeightedPoolTradingCurve, alpha, beta; pathological)
    if beta == Inf
        return Inf
    end
    (beta, dydp_p(tc, beta) * beta)
end

