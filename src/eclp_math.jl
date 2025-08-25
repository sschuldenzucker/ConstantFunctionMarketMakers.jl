"""
Module for the ECLP math. This is kept separate from the superstructure in `ConstantFunctionMarketMakers` to
avoid namespace pollution.
"""
module ECLPMath

using Test
using LinearAlgebra
using NLsolve
using Symbolics
using Nemo  # required for symbolic solve

"""
Master params structure. Note that this stores a bit too many params for some operations, e.g.,
alpha and beta are not required for tau etc.

SOMEDAY We may want a separate cache struct.
"""
struct Params
    alpha::Float64
    beta::Float64
    q::Float64
    lam::Float64

    # Implied
    s::Float64
    c::Float64

    function Params(alpha, beta, q, lam)
        @assert 0 < alpha < q < beta < Inf
        @assert 1 <= lam
        s, c = sc_from_q(q)
        new(alpha, beta, q, lam, s, c)
    end
end

function A_mat(m::Params)
    c, s, lam = m.c, m.s, m.lam
    [
        c/lam -s/lam
        s c
    ]
end

function A_inv_mat(m::Params)
    c, s, lam = m.c, m.s, m.lam
    [
        c*lam s
        -s*lam c
    ]
end

# TODO bunch of sanity asserts, e.g., prices within bounds.

# TODO are we gonna use SVector / SPair or just arrays? Would have to decide at some point I think.
eta(pc) = 1 / sqrt(1 + pc^2) * [pc, one(pc)]

sc_from_q(q) = eta(q)

function zeta(m::Params, p)
    ap = A_mat(m) * [-1, p]
    dot([0, 1], ap) / (-dot([1, 0], ap))
end

tau(m::Params, p) = eta(zeta(m, p))

# Sanity check
# NB alpha and beta are arbitrary to pass checks
@test tau(Params(0.5, 1.5, 1, 2), 1) == [0, 1]

A_inv_tau(m::Params, p) = A_inv_mat(m) * tau(m, p)

function chi(m::Params)
    # SOMEDAY A_inv is cached here, but tau recomputes A_mat, which is sad.
    [A_inv_tau(m, m.beta)[1], A_inv_tau(m, m.alpha)[2]]
end

"""
Residual from the ECLP invariant equation.  
"""
function residual_eclp(m::Params, t, l)
    norm_sq = t1 -> dot(t1, t1)
    # SOMEDAY another instance of recalculation of A_mat.
    # Maybe cache in params.
    norm_sq(A_mat(m) * (t - l * chi(m))) - l^2
end

function t_l_p(m::Params, l, p)
    l * (chi(m) - A_inv_tau(m, p))
end

"""Price as a fct of r=r(x, y), x, y. See the paper and calcSpotPrice0in1() in GyroECLPMath.sol."""
function p_l_t(m::Params, l, t)
    # I'm honestly not quite sure what these incantations do but it worked fine in the Sage version
    # so I guess it's fine.
    A = A_mat(m)
    tpp = A * (t - l * chi(m))
    pc = [tpp[1] / tpp[2], 1]
    pgx = pc' * (A * [1, 0])
    px = pgx / (pc' * (A * [0, 1]))
    px
end

function t_plus(m::Params, l)
    # SOMEDAY could be cached in Params if we need it a lot.
    ma = [
        1.0 0
        0.0 -1.0
    ]
    xyplus_1 = ma * A_inv_mat(m) * (tau(m, m.beta) - tau(m, m.alpha))
    l .* xyplus_1
end

"""
Calculate the invariant numerically  
"""
# TODO run a simple test case to check this actually is fine.
function find_invariant(m::Params, t)
    res = nlsolve(v -> residual_eclp(m, t, v[1]), [1.0]; autodiff = :forward)
    @assert res.f_converged
    res.zero[1]
end

# TODO move these to some variables module. They are used in stableswap too.
_x = Symbolics.variable("x")
_y = Symbolics.variable("y")
_p = Symbolics.variable("p")

"""
Returns a lambda that computes x(y). Uses symbolics under the hood (but not exposed).

The lambda takes one argument (namely, y) and returns a Float64.
"""
function mk_x_y(m::Params, l)
    slns = symbolic_solve(residual_eclp(m, [_x, _y], l), _x)
    sln = minimum(slns)
    fn = eval(build_function(sln, _y))
    y -> Float64(Base.invokelatest(fn, y))
end

"""
See `mk_x_y`
"""
function mk_y_x(m::Params, l)
    slns = symbolic_solve(residual_eclp(m, [_x, _y], l), _y)
    sln = minimum(slns)
    fn = eval(build_function(sln, _x))
    x -> Float64(Base.invokelatest(fn, x))
end

"""
Liquidity density. Otherwise like `mk_x_y`.
"""
function mk_dydp_p(m::Params, l)
    y_p = t_l_p(m, l, _p)[2]
    expr = expand_derivatives(Differential(_p)(y_p))
    fn = eval(build_function(expr, _p))
    p -> Float64(Base.invokelatest(fn, p))
end

end
