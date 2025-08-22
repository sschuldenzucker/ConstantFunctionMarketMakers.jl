"""
Module for the StableSwap math. This is kept separate from the superstructure in `CalibrationTools` to
avoid namespace pollution.

Note that, in contrast to most other AMMs, we don't have formulas based on p as the input variable, only on y.

Note that the CLP only has one parameter, the amplification factor A. We take the "math form" of the StableSwap. There is also the "code form", where the meaning of the A parameter differs by a factor (2 for n=2).

SOMEDAY fix that, it's confusing.

In contrast to the SageMath version of this code (stableswap.sage), this file only implements stableswap for 2 assets. This restriction can still be solved analytically.
"""
# Code structure is parallel to `eclp_math.jl` as much as possible.
module StableSwapMath

using Test
using LinearAlgebra
using NLsolve
using Symbolics, Nemo
using ForwardDiff

using ..Variables

"""
Residual from the StableSwap invariant equation.

Works for any number of assets, but has to be indicated explicitly for now. (SOMEDAY just use n=length(t))
"""
function residual_stableswap(a, t, l; n = 2)
    a * n^n * sum(t) + l - (a * l * n^n + l^(n + 1) / (n^n * prod(t)))
end

function find_invariant(a, t; n = 2)
    l0 = sum(t)  # sum(t) is the correct invariant if all of t's elements are the same.
    res = nlsolve(v -> residual_stableswap(a, t, v[1]; n), [l0]; autodiff = :forward)
    @assert res.f_converged
    res.zero[1]
end

# SOMEDAY Can we have `p_l_t` analytically (like, an actual formula)?

# We symbolically solve the stableswap equation (n=2) for _general_ variables (without replacing a
# first; fully static)
# We then operate no the resulting compiled functions.

# NB we don't have x_y_d etc. right now but we could have them easily.

const _y_x_gen = let
    slns = symbolic_solve(expand(residual_stableswap(_a, [_x, _y], _l) * _x * _y), _y)
    # NB I just know that slns[1] is the valid one (and the other one is always negative, thus invalid.)
    # TODO whether it's slns[1] or slns[2] is a bit flaky. But I don't wanna use max here to keep this differentiable.
    sln = simplify(slns[2])
    eval(build_function(sln, _a, _l, _x))
end

const _x_y_gen = let
    slns = symbolic_solve(expand(residual_stableswap(_a, [_x, _y], _l) * _x * _y), _x)
    # NB I just know that slns[1] is the valid one (and the other one is always negative, thus invalid.)
    # TODO whether it's slns[1] or slns[2] is a bit flaky. But I don't wanna use max here to keep this differentiable. Ideally we could add a constraint on x.
    sln = simplify(slns[2])
    eval(build_function(sln, _a, _l, _y))
end

# SOMEDAY should we typecast to Float64 (from BigFloat) here?
function y_x_gen(a, l, x)
    if x == 0.0
        return Inf
    end
    if x == Inf
        return 0.0
    end
    _y_x_gen(a, l, x)
end
function x_y_gen(a, l, y)
    if y == 0.0
        return Inf
    end
    if y == Inf
        return 0.0
    end
    _x_y_gen(a, l, y)
end

# The price at y is the derivative of y(x) at x(y).
function p_y_gen(a, l, y)
    # Special cases, not sure these are gonna be useful.
    if y == 0.0
        return 0.0
    end
    if y == Inf
        return Inf
    end

    x = x_y_gen(a, l, y)
    @debug "p_y_gen evaluation" y x
    -ForwardDiff.derivative(x -> y_x_gen(a, l, x), x)
end

# dydp(y) = 1 / dpdy(y)
# SOMEDAY is this kinda slow?
dydp_y_gen(a, l, y) = 1.0 / ForwardDiff.derivative(y -> p_y_gen(a, l, y), y)

# p_x_gen(a, l, x) = -ForwardDiff.derivative(x -> y_x_gen(a, l, x), x)

@test y_x_gen(10.0, 200.0, 100.0) == 100.0
@test p_y_gen(10.0, 200.0, 100.0) == 1.0

function find_t_l_p(a, l, p)
    # Special cases. NB the numerics of this are gonna be unstable way before this point, but it's
    # useful being able to put these for use with `ConcentratedAMM`.
    if p == 0.0
        return [Inf, 0.0]
    end
    if p == Inf
        return [0.0, Inf]
    end

    y0 = l / 2  # exact solution for p=1
    @debug "find_t_l_p" y0

    # Transform the input so that we don't go to non-positive values.
    # SOMEDAY we should prob use a different solver here instead. This transformation distorts the
    # error a bit. Could lead to (bogus non-convergence) or (excessive error for very small y).
    # Maybe we should use the relative error in the inner lambda to compensate?
    trans(x) = exp(x)
    untrans(x) = log(x)
    y0u = untrans(y0)
    res = nlsolve(v -> p_y_gen(a, l, trans(v[1])) - p, [y0u]; autodiff = :forward)
    @assert res.f_converged
    y = trans(res.zero[1])

    [x_y_gen(a, l, y), y] .|> Float64
end

# -------------------------------------
# DEPRECATED and not used.
# It's usable but I've found the `_gen` variants above to be better. Probably.
# -------------------------------------

"""
Struct that stores a bunch of computed (and cached for speed) lambdas that come from analytic solutions (for n=2).
"""
# NB having a bunch of lambdas here is normally a code smell but there doesn't seem any useful data
# to base this on. Maybe in the future use the formuals and then maybe we can get rid of some type
# issues but idk.
struct Cache
    x_y::Function
    y_x::Function
    p_y::Function
    dydp_y::Function
    # We do *not* have dydp_p! (at least not fast)
end

function Cache(a, l; n = 2)
    @assert n == 2 "Not implemented for n > 2"

    # NOTE simplify() is really important here. O/w error accumulation sends you to hell.
    # I got complete *nonsense* for p_y without it.
    # SOMEDAY maybe I should use AD to get p_y instead of the symbolic derivatives.
    # SOMEDAY is simplify() really slow tough?

    x_y_fml = let
        # manual expansion to get a polynomial as of https://github.com/JuliaSymbolics/Symbolics.jl/issues/1268
        slns = symbolic_solve(expand(residual_stableswap(a, [_x, _y], l) * _x * _y), _x)
        # NB I just know that slns[1] is the valid one (and the other one is always negative, thus invalid.)
        simplify(slns[1])
    end
    @info "x_y_fml"
    y_x_fml = let
        slns = symbolic_solve(expand(residual_stableswap(a, [_x, _y], l) * _x * _y), _y)
        simplify(slns[1])
    end
    @info "y_x_fml"
    # p_y_fml = 1 / expand_derivatives(Differential(_y)(x_y_fml))
    p_y_fml = simplify(
        substitute(-expand_derivatives(Differential(_x)(y_x_fml)), Dict(_x => x_y_fml)),
    )
    @info "p_y_fml"

    # simplify THIS is the slow one!! (not doing that anymore; not sure if that's bad)
    # dydp_y_fml = simplify(1 / expand_derivatives(Differential(_y)(p_y_fml)))
    # dydp_y_fml = 1 / simplify(expand_derivatives(Differential(_y)(p_y_fml)))
    dydp_y_fml = 1 / expand_derivatives(Differential(_y)(p_y_fml))
    @info "dydp_y_fml"

    # TODO this leads to a world age issue that prevents specialization and requires invokelatest below.
    # Instead, eval it into a specific module at global scope.
    # I think my version of Symbolics also defines an *anonymous* function, unclear wat do.
    # This is relevant here b/c of automatic differentiation in `find_t_l_p`.
    # See https://chatgpt.com/share/6804ebf2-93a4-800f-a691-73498b1d34a9 (I think this is a bit funky, maybe wrong)
    x_y_fn = eval(build_function(x_y_fml, _y))
    y_x_fn = eval(build_function(y_x_fml, _x))
    p_y_fn = eval(build_function(p_y_fml, _y))
    dydp_y_fn = eval(build_function(dydp_y_fml, _y))

    @info "Built functions"

    Cache(
        # y -> Float64(x_y_fn(y)),
        # x -> Float64(y_x_fn(x)),
        # y -> Float64(p_y_fn(y)),
        # y -> Float64(dydp_y_fn(y)),
        y -> Base.invokelatest(x_y_fn, y),
        x -> Base.invokelatest(y_x_fn, x),
        y -> Base.invokelatest(p_y_fn, y),
        y -> Base.invokelatest(dydp_y_fn, y),
    )
end

function find_t_l_p(cache::Cache, l, p)
    y0 = l / 2  # exact solution for p=1
    @info cache.p_y(y0)
    res = nlsolve(v -> cache.p_y(v[1]) - p, [y0]; autodiff = :forward)
    @info res
    @assert res.f_converged
    y = res.zero[1]

    [cache.x_y(y), y] .|> Float64
end

end
