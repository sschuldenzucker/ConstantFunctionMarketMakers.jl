module CalibrationTools

export AMM,
    kind_name,
    short_str,
    flip,
    t_l_p,
    p_l_t,
    residual,
    l_t,
    mk_tc,
    tc_t,
    tc_v_p,
    tc_l_p,
    trading_curve_from_balances,
    trading_curve_from_value_price,
    trading_curve_from_invariant_price,
    t_l_p,
    p_l_t,
    alphabeta,
    p_y,
    dydp_y,
    sample_dydp_p,
    sample_dydlogp_p,
    t_plus,
    v_init,
    dydp_max,
    dydlogp_max

using LinearAlgebra
# SOMEDAY for even more performance (but prob doesn't matter)
# using StaticArrays
using Test

using Accessors

using ForwardDiff

using ParametricAdaptiveSampling

include("util.jl")
using .Util

# Lets modules do `using ..Variables` to bring _x etc. into scope.
include("variables.jl")

# High-level infra and method definitions for AMMs
include("infra.jl")

# AMMs
include("cpmm.jl")
include("c2lp.jl")
include("eclp.jl")
include("stableswap.jl")
include("weighted.jl")
include("cpsmm.jl")

# TODO constant ld and constant relative LD
# TODO flip operation (view from the other side) - Should this be an AMM wrapper or a change in parameters? I guess either, depending on the AMM.

# External (exported) tools
include("tools.jl")

# AMM extensions
include("rate_scaled_amm.jl")
include("amm_with_meta.jl")

# Algorithms
include("concentrate.jl")

# TODO and all of the other tools as well (sampling, plotting after sampling, handling StableSwap vs ECLP, handling bounds suggestions)
#
# For plotting and sampling, I'll need a trait specifying how to get the LD (dydp_p or dydp_y for StableSwap). This can default to hasmethod(dydp_p) but needs to be hand-written for the RateScaledAMM. Also ill-defined for (e.g.) sums of curves. Maybe we want a higher-level structure.

# TODO more types:
# - TradingCurve as a sum of TCs. Maybe CFMM as a sum of CFMMs (I guess). The use case here is combining LD curves, parallel pools, and (as a special case of that) Uni v3 TCs (TradingCurve only, CFMM not defined).
# - Trading through a chain of pools. (this is a bit unclear)

# --------------------------------------
# Late Infra
# --------------------------------------

end
