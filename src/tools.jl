
export union_priceset, intersect_priceset

"""
    union_priceset(amms) :: Vector{Tuple{Float64,Float64}}

Return the union of the price sets of a set of AMMs as a minimal union of intervals. Can be used with both AMMs and TCs.
"""
union_priceset(amms) = union_intervals(alphabeta.(amms))
union_priceset(amm::Union{AMM,TradingCurve}) = [alphabeta(amm)]

"""
    intersect_priceset(amms) :: Union{Nothing,Tuple{Float64,Float64}}

Return the intersection of the price sets of a set of AMMs. Returns Nothing if the AMMs don't intersect. Can be used with both AMMs and TCs.
"""
intersect_priceset(amms) = intersect_intervals(alphabeta.(amms))
intersect_priceset(amm::Union{AMM,TradingCurve}) = [alphabeta(amm)]

function union_intervals(ivls::Vector{Tuple{T,S}}) where {T<:Real,S<:Real}
    if isempty(ivls)
        return []
    end

    ivls = sort(ivls; by = first)

    cur = first(ivls)
    ret = []
    for ivl in Iterators.rest(ivls)
        (a1, b1) = ivl
        (a, b) = cur
        if a <= a1 <= b
            cur = (a, max(b1, b))
        else
            push!(ret, cur)
            cur = ivl
        end
    end
    push!(ret, cur)
    ret
end

@test union_intervals([(0, 1), (2, 4), (0.5, 3)]) == [(0, 4)]
@test union_intervals([(0, 1), (4, 5), (0.5, 3), (4.5, 6)]) == [(0, 3), (4, 6)]
@test union_intervals([(2, Inf), (1, 2), (0, 0.5)]) == [(0, 0.5), (1, Inf)]

function intersect_intervals(ivls::Vector{Tuple{T,S}}) where {T<:Real,S<:Real}
    if isempty(ivls)
        return nothing
    end

    a = maximum(first, ivls)
    b = minimum(last, ivls)

    if a <= b
        (a, b)
    else
        nothing
    end
end

@test intersect_intervals([(0.5, 2), (0, 1), (0.7, Inf)]) == (0.7, 1)
@test intersect_intervals([(0.5, 2), (2, 3)]) == (2, 2)
@test intersect_intervals([(0.5, 2), (2.5, 3)]) |> isnothing

