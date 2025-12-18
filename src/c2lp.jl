export C2LP, C2LPTradingCurve

module C2LPMath

struct Params
    alpha::Float64
    beta::Float64
    function Params(alpha, beta)
        # SOMEDAY technically having them infinite is probably fine? (though our impl may break)
        @assert 0 < alpha < beta < Inf
        new(alpha, beta)
    end
end

calc_ab(m::Params, l = 1.0) = [l / sqrt(m.beta), l * sqrt(m.alpha)]

function residual_c2lp(m::Params, t, l)
    ab = calc_ab(m, l)
    prod(t .+ ab) - l^2
end

function find_invariant(m::Params, t)
    res = nlsolve(v -> residual_c2lp(m, t, v[1]), [1.0]; autodiff = :forward)
    @assert res.f_converged
    return res.zero[1]
end

t_plus(m::Params, l = 1.0) =
    l .* [1.0 / sqrt(m.alpha) - 1.0 / sqrt(m.beta), sqrt(m.beta) - sqrt(m.alpha)]

t_l_p(m::Params, l, p) =
    [l * (1.0 / sqrt(p) - 1.0 / sqrt(m.beta)), l * (sqrt(p) - sqrt(m.alpha))]

function p_l_t(m::Params, l, t)
    ab = calc_ab(m, l)
    tv = t .+ ab
    tv[2] / tv[1]
end

function x_y(m::Params, l, y)
    a, b = calc_ab(m, l)
    l^2 / (y + b) - a
end

function y_x(m::Params, l, x)
    a, b = calc_ab(m, l)
    l^2 / (x + a) - b
end

end  # module C2LPMath

import .C2LPMath

@doc raw"""
    C2LP(alpha, beta)

Gyroscope 2CLP, aka concentrated constant-product, aka a single Uniswap v3 bracket.

```math
(x + a) \cdot (y + b) = L^2
```

where $a = l / \sqrt{\beta}$ and $b = l \cdot \sqrt{\alpha}$.

See also the [Gyroscope docs](https://docs.gyro.finance/pools/2-clps.html).
"""
struct C2LP <: AMM
    params::C2LPMath.Params
end

C2LP(alpha, beta) = C2LP(C2LPMath.Params(alpha, beta))

flip(amm::C2LP) = C2LP(C2LPMath.Params(1.0 / amm.params.beta, 1.0 / amm.params.alpha))

function display_params(amm::C2LP; as_bp::Bool)
    # SOMEDAY should we just make this choice automatic?
    # SOMEDAY other display styles, e.g., FDV (relative to something indicated in a union type)?
    mogrify = as_bp ? fmt_todev : fmt_auto
    ["α" => amm.params.alpha |> mogrify, "β" => amm.params.beta |> mogrify]
end

kind_name(::C2LP) = "2CLP"

alphabeta(amm::C2LP) = (amm.params.alpha, amm.params.beta)

residual(amm::C2LP, l, t) = C2LPMath.residual_c2lp(amm.params, t, l)

l_t(amm::C2LP, t) = C2LPMath.find_invariant(amm.params, t)

t_l_p(amm::C2LP, l, p) = C2LPMath.t_l_p(amm.params, l, p)
p_l_t(amm::C2LP, l, t) = C2LPMath.p_l_t(amm.params, l, t)

mk_tc(common::TCCommon{C2LP}) = C2LPTradingCurve(common)

struct C2LPTradingCurve <: TradingCurve
    common::TCCommon{C2LP}
end

x_y(tc::C2LPTradingCurve, y) = C2LPMath.x_y(tc.common.amm.params, tc.common.l, y)
y_x(tc::C2LPTradingCurve, x) = C2LPMath.y_x(tc.common.amm.params, tc.common.l, x)

t_plus(tc::C2LPTradingCurve) = C2LPMath.t_plus(tc.common.amm.params, tc.common.l)

dydp_p(tc::C2LPTradingCurve, p) = ForwardDiff.derivative(p -> t_p(tc, p)[2], p)

_dydp_max(tc::C2LPTradingCurve, alpha, beta) = (alpha, dydp_p(tc, alpha))
_dydlogp_max(tc::C2LPTradingCurve, alpha, beta) =
    let p = beta
        (p, dydp_p(tc, p) * p)
    end

