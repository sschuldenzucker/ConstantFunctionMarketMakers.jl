
export SumTradingCurve, sum_tcs

struct SumTradingCurve <: TradingCurve
    tcs::Vector{TradingCurve}
end

sum_tcs(tcs) = SumTradingCurve(tcs)

kind_name(::SumTradingCurve) = "Sum"

short_str(tc::SumTradingCurve; as_bp::Bool) = join(short_str.(tc.tcs; as_bp), " + ")

# TODO not quite right: there can be holes! Needs some work.
alphabeta(tc::SumTradingCurve) = envelop_priceset(tc.tcs)

flip(tc::SumTradingCurve) = SumTradingCurve(flip.(tc.tcs))

# Component-wise sum
t_plus(tc::SumTradingCurve) = foldl((a, b) -> a .+ b, alphabeta.(tc.tcs); init = (0, 0))

# NOTE: These are gonna be slow for StableSwap but there's not much I can do: summing _fundamentally_ works across prices!

dydp_p(tc::SumTradingCurve, p) =
    sum(tc.tcs) do td
        is_in_range(td, p) ? dydp_p(td, p) : 0.0
    end

t_p(tc::SumTradingCurve, p) = sum(td -> t_p(td, p), tc.tcs)

"""
    _make_bracket_nondegenerate(f, a, b)

From a potentially infinite (from 0.0 OR to Inf) bracket for f, create a finite one.
"""
function _make_bracket_nondegenerate(f, a, b; stepsize = 2.0)
    if a == 0.0
        a = b < Inf ? b / stepsize : 1.0
        while f(a) >= 0
            a /= stepsize
        end
    end

    if b == Inf
        b = a * stepsize
        while f(b) <= 0
            b *= stepsize
        end
    end

    a, b
end

using Roots

# Note that we cannot just use the `p_x` of `tc.tcs` even when they are fast. This is not how it works.
function p_x(tc::SumTradingCurve, x)
    alpha, beta = alphabeta(tc)

    if x == 0.0
        return beta
    end
    if x == Inf
        return alpha
    end

    # This function is monotonic and has a root.
    f(p) = x - x_p(tc, p)

    a, b = _make_bracket_nondegenerate(f, alpha, beta)
    find_root(f, (a, b), Roots.A42())
end

function p_y(tc::SumTradingCurve, y)
    alpha, beta = alphabeta(tc)

    if y == 0.0
        return alpha
    end
    if y == Inf
        return beta
    end

    # This function is monotonic and has a root.
    f(p) = y_p(tc, p) - y

    a, b = _make_bracket_nondegenerate(f, alpha, beta)
    find_root(f, (a, b), Roots.A42())
end

x_y(tc::SumTradingCurve, y) = x_p(tc, p_y(tc, y))
y_x(tc::SumTradingCurve, x) = y_p(tc, p_x(tc, x))

"""
For `SumTradingCurve`, we sample parametrically along the price dimension b/c that's more efficient than sampling `y_x`, which is quite expensive.
"""
function sample_x_y(tc::SumTradingCurve, xmin, xmax = x_plus(tc))
    pmin = p_x(tc, xmin)
    pmax = p_x(tc, xmax)
    sample_adaptive_parametric(t_p, pmin, pmax)
end

# We do NOT overload _dydp_max and _dydlogp_max b/c we can't do anything clever here. This can be slow, though, b/c it runs a bunch of numerical searches within each other.

