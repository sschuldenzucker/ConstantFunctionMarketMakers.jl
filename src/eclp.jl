
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

l_t(eclp::ECLP, t) = ECLPMath.find_invariant(eclp.params, t)

t_l_p(eclp::ECLP, l, p) = ECLPMath.t_l_p(eclp.params, l, p)
p_l_t(eclp::ECLP, l, t) = ECLPMath.p_l_t(eclp.params, l, t)

function mk_tc(common::TCCommon{ECLP})
    params = common.amm.params
    ECLPTradingCurve(
        common,
        # SOMEDAY we could replace these with the generic solution method used for StableSwap. Should be faster & more stable.
        ECLPMath.mk_x_y(params, common.l),
        ECLPMath.mk_y_x(params, common.l),
        ECLPMath.mk_dydp_p(params, common.l),
        ECLPMath.t_plus(params, common.l),
    )
end

struct ECLPTradingCurve <: TradingCurve
    common::TCCommon{ECLP}

    # Cached results
    _x_y::Function
    _y_x::Function
    _dydp_p::Function
    _t_plus::Vector{Float64}
end

x_y(tc::ECLPTradingCurve, y) = tc._x_y(y)
y_x(tc::ECLPTradingCurve, x) = tc._y_x(x)
dydp_p(tc::ECLPTradingCurve, p) = tc._dydp_p(p)

t_p(tc::ECLPTradingCurve, p) = t_l_p(tc.common.amm, tc.common.l, p)

t_plus(tc::ECLPTradingCurve) = tc._t_plus

