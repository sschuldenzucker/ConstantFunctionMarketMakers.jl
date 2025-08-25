
include("./eclp_math.jl")
import .ECLPMath

export ECLP, ECLPTradingCurve

struct ECLP <: AMM
    params::ECLPMath.Params
end

ECLP(alpha, beta, q, lam) = ECLP(ECLPMath.Params(alpha, beta, q, lam))

flip(amm::ECLP) = ECLP(
    ECLPMath.Params(
        1.0 / amm.params.beta,
        1.0 / amm.params.alpha,
        1.0 / amm.params.q,
        amm.params.lam,
    ),
)

function display_params(eclp::ECLP; as_bp::Bool)
    # SOMEDAY should we just make this choice automatic?
    # SOMEDAY other display styles, e.g., FDV (relative to something indicated in a union type)?
    mogrify = as_bp ? fmt_todev : fmt_auto
    [
        "α" => eclp.params.alpha |> mogrify,
        "β" => eclp.params.beta |> mogrify,
        "peg" => eclp.params.q |> mogrify,
        "λ" => eclp.params.lam |> fmt_auto,
    ]
end

alphabeta(eclp::ECLP) = (eclp.params.alpha, eclp.params.beta)

residual(eclp::ECLP, l, t) = ECLPMath.residual_eclp(eclp.params, t, l)

l_t(eclp::ECLP, t) = ECLPMath.calc_invariant(eclp.params, t)

t_l_p(eclp::ECLP, l, p) = ECLPMath.t_l_p(eclp.params, l, p)
p_l_t(eclp::ECLP, l, t) = ECLPMath.p_l_t(eclp.params, l, t)

function mk_tc(common::TCCommon{ECLP})
    params = common.amm.params
    ECLPTradingCurve(common, ECLPMath.t_plus(params, common.l))
end

struct ECLPTradingCurve <: TradingCurve
    common::TCCommon{ECLP}

    # Cached results
    _t_plus::Vector{Float64}
end

x_y(tc::ECLPTradingCurve, y) = ECLPMath.x_y(tc.common.amm.params, tc.common.l, y)
y_x(tc::ECLPTradingCurve, x) = ECLPMath.y_x(tc.common.amm.params, tc.common.l, x)
dydp_p(tc::ECLPTradingCurve, p) = ForwardDiff.derivative(p -> t_p(tc, p)[2], p)  # The generic implementation

t_plus(tc::ECLPTradingCurve) = tc._t_plus

