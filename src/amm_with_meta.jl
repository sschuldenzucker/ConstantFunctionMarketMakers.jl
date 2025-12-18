
export PlotMeta, AMMWithMeta, TradingCurveWithMeta, named

struct PlotMeta
    name::Union{Nothing,String}
    plot_attr::Dict{Symbol,Any}
end

PlotMeta(name) = PlotMeta(name, Dict())

"""
An AMM with metadata of type `M` attached. Often, `M = PlotMeta`.
"""
struct AMMWithMeta{M,A<:AMM} <: AMM
    amm::A
    meta::M
end

flip(amm::AMMWithMeta) = AMMWithMeta(flip(amm.amm), amm.meta)
flip(amm::AMMWithMeta{PlotMeta}) =
    AMMWithMeta(flip(amm.amm), PlotMeta(amm.meta.name * " (flipped)", amm.meta.plot_attr))

kind_name(amm::AMMWithMeta) = kind_name(amm.amm)

function short_str(amm::AMMWithMeta{PlotMeta,A}; as_bp::Bool) where {A}
    ret = short_str(amm.amm; as_bp)
    if !isnothing(amm.meta.name)
        ret = amm.meta.name * " " * ret
    end
    ret
end

alphabeta(amm::AMMWithMeta) = alphabeta(amm.amm)

residual(amm::AMMWithMeta, l, t) = residual(amm.amm, l, t)

l_t(amm::AMMWithMeta, l) = l_t(amm.amm, l)

t_l_p(amm::AMMWithMeta, l, p) = t_l_p(amm.amm, l, p)

p_l_t(amm::AMMWithMeta, l, p) = p_l_t(amm.amm, l, p)

function mk_tc(common::TCCommon{AMMWithMeta{M,A}}) where {M,A}
    common_inner = TCCommon(common.amm.amm, common.l, common.t_init, common.p_init)
    TradingCurveWithMeta(common, mk_tc(common_inner))
end

"""
A TradingCurve with metadata of type `M` attached. Often, `M = PlotMeta`.
"""
struct TradingCurveWithMeta{M,A<:AMM,T<:TradingCurve} <: TradingCurve
    common::TCCommon{AMMWithMeta{M,A}}
    inner::T
end

t_plus(tc::TradingCurveWithMeta) = t_plus(tc.inner)

x_y(tc::TradingCurveWithMeta, y) = x_y(tc.inner, y)
y_x(tc::TradingCurveWithMeta, x) = y_x(tc.inner, x)
dydp_p(tc::TradingCurveWithMeta, p) = dydp_p(tc.inner, p)
dydp_y(tc::TradingCurveWithMeta, y) = dydp_y(tc.inner, y)
_sample_dydp_p(tc::TradingCurveWithMeta, alpha, beta) =
    _sample_dydp_p(tc.inner, alpha, beta)
_sample_dydlogp_p(tc::TradingCurveWithMeta, alpha, beta) =
    _sample_dydlogp_p(tc.inner, alpha, beta)
_dydp_max(tc::TradingCurveWithMeta, alpha, beta; pathological) =
    _dydp_max(tc.inner, alpha, beta; pathological)
_dydlogp_max(tc::TradingCurveWithMeta, alpha, beta; pathological) =
    _dydlogp_max(tc.inner, alpha, beta; pathological)

# Convenience method
"""
    named(amm, name) <: AMM
    amm |> named(name)

Attach a custom name to an AMM
"""
named(amm::AMM, name::String) = AMMWithMeta(amm, PlotMeta(name))
named(amm::AMMWithMeta{PlotMeta,A}, name) where {A} = @set amm.meta.name = name
named(name::String) = amm -> named(amm, name)

# TODO as_bp is really ill-defined. This should all be automatic, see util.jl.
# and as_bp shouldn't be there.
plot_attr(amm::AMMWithMeta{PlotMeta,A}) where {A} =
    merge(amm.meta.plot_attr, Dict(:label => short_str(amm; as_bp = false)))
plot_attr(tc::TradingCurveWithMeta{PlotMeta,A,T}) where {A,T} = plot_attr(tc.common.amm)

