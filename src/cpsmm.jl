
export CPSMM, CPSMMTradingCurve

"""
A constant-power-sum-pool (like YieldSpace) following \$x^delta + y^delta = L^delta\$

where delta is strictly between 0 and 1. delta -> 1 is the constant-sum AMM and delta -> 0 is the CPMM.
"""
struct CPSMM <: AMM
    delta::Float64

    function CPSMM(delta)
        @assert 0 < delta < 1
        new(delta)
    end
end

flip(amm::CPSMM) = amm

display_params(cpsmm::CPSMM; as_bp::Bool) = ["δ" => fmt_auto(cpsmm.delta)]

alphabeta(::CPSMM) = (0.0, Inf)

residual(cpsmm::CPSMM, l, t) = sum(t .^ cpsmm.delta) - l^cpsmm.delta

l_t(cpsmm::CPSMM, t) = (t[1]^cpsmm.delta + t[2]^cpsmm.delta)^(1 / cpsmm.delta)

function t_l_p(cpsmm::CPSMM, l, p)
    delta = cpsmm.delta
    theta = (1.0 / (1.0 - delta) - 1.0)
    l .* [(1.0 + p^theta)^(-1.0 / delta), (1.0 + p^(-theta))^(-1.0 / delta)]
end

p_l_t(cpsmm::CPSMM, l, t) = ((t[2] / t[1])^(1.0 - cpsmm.delta))

mk_tc(common::TCCommon{CPSMM}) = CPSMMTradingCurve(common)

struct CPSMMTradingCurve <: TradingCurve
    common::TCCommon{CPSMM}
end

function x_y(tc::CPSMMTradingCurve, y)
    delta = tc.common.amm.delta
    (tc.common.l^delta - y^delta)^(1.0 / delta)
end

function y_x(tc::CPSMMTradingCurve, x)
    delta = tc.common.amm.delta
    (tc.common.l^delta - x^delta)^(1.0 / delta)
end

t_plus(tc::CPSMMTradingCurve) =
    let l = tc.common.l
        (l, l)
    end

# SOMEDAY we actually know this analytically, see the paper. (but ForwardDiff should honestly be fine, too)
dydp_p(tc::CPSMMTradingCurve, p) = ForwardDiff.derivative(p -> t_p(tc, p)[2], p)

using Optim

function _dydp_max(tc::CPSMMTradingCurve, alpha, beta; pathological)
    # We know that the maximum is gonna be <= 1.0, and we give it a bit of slack
    if beta > 2.0
        beta = max(2.0 * alpha, 2.0)
    end

    # NB this bottoms out at around eps() ~ 1e-16. We could try to work around the machine precision but it doesn't seem to be worth it.
    # This works for alpha = 0.0. We know that the maximum is gonna be <= 1.0
    res = optimize(p -> -dydp_p(tc, p), alpha, beta)
    @assert res.stopped_by.converged
    (res.minimizer, -res.minimum)
end

function _dydlogp_max(tc::CPSMMTradingCurve, alpha, beta; pathological)
    function go(pmin, pmax)
        res = optimize(p -> -dydp_p(tc, p) * p, pmin, pmax)
        # res = optimize(p -> -cpsmm_ld_analytic(delta, p) * p, pmin, pmax)
        @assert res.stopped_by.converged
        (res.minimizer, -res.minimum)
    end

    # We know that the maximum is going to be >= 1.0 and give it a bit of slack.
    if alpha < 0.5
        alpha = min(0.5 * beta, 0.5)
    end

    # This avoids some numerical issues I don't fully understand.
    if tc.common.amm.delta >= 0.5
        # I know that the maximum is going to be well below 5.
        beta = max(2.0 * alpha, 5.0)
        return go(alpha, beta)
    end

    # We need to make beta finite and 1e100 is approximately the point where numerics become a real problem.
    # NB This is gonna clip at 1e100 but we tolerate that, doesn't seem worth fixing.
    # I tried unconstrained optimization with some transforms but wasn't successful for some reason.
    beta = min(beta, 1e100)
    go(alpha, beta)
end
