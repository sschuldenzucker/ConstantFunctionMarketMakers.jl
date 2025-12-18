export CPMM, CPMMTradingCurve

@doc raw"""
    CPMM()

Constant-product AMM.

```math
x \cdot y = L^2
```
"""
struct CPMM <: AMM end

flip(amm::CPMM) = amm

residual(::CPMM, l, t) = prod(t) - l^2

l_t(::CPMM, t) = sqrt(prod(t))

t_l_p(::CPMM, l, p) = [l / sqrt(p), l * sqrt(p)]

p_l_t(::CPMM, l, t) =
    let (x, y) = t
        y / x
    end

mk_tc(common::TCCommon{CPMM}) = CPMMTradingCurve(common)

struct CPMMTradingCurve <: TradingCurve
    common::TCCommon{CPMM}
end

t_plus(::CPMMTradingCurve) = (Inf, Inf)

x_y(tc::CPMMTradingCurve, y) = tc.common.l^2 / y
y_x(tc::CPMMTradingCurve, x) = tc.common.l^2 / x

dydp_p(tc::CPMMTradingCurve, p) = 0.5 * tc.common.l / sqrt(p)

# NB this works just fine with 0 and inf
_dydp_max(tc::CPMMTradingCurve, alpha, beta; pathological) = (alpha, dydp_p(tc, alpha))
_dydlogp_max(tc::CPMMTradingCurve, alpha, beta; pathological) =
    let p = beta
        (p, dydp_p(tc, p) * p)
    end

