module Util

export todev, fromdev, fmt_todev, fmt_auto, intersect_intervals, norm_sq

using Printf
using LinearAlgebra

todev(x) = (x - one(x)) * 1e4
fromdev(x) = one(x) + x * 1e-4

"Format as dev, e.g., 0.9990 becomes '-10bp'"
fmt_todev(x) =
    let xd = todev(x)
        @sprintf("%+.0f", xd) * "bp"
    end

"""
    fmt_auto(x; digits = 2, suffixes = true, bp = true, bp1 = true)

Format using an automatic number of decimals based on the size of the number.

digits: number of digits of precision to show
suffixes: if true, use suffixes 'M', 'k' for large numbers where appropriate.
bp: if true, write small numbers as bp
dev: if true, write numbers close to 1 as +/- bp
"""
function fmt_auto(x; digits = 2, suffixes = true, bp = true, dev = true)
    # TODO not implemented, got bored
    # TODO doesn't support negative x right now. Split into abs and sign and multiply later.
    # TODO also a bit wonky I think
    # Here we always show 4 digits.
    # if dev && abs(x - one(x)) < 100
    #     return fmt_todev(x)
    # end

    # if bp && x < 1e-2
    #     return @sprintf("%.0fbp", x * 1e4)
    # end

    # oom = log10(x) |> floor |> Int
    # digits_factor = 10.0^digits

    # if suffixes && x >= 1e6 * digits_factor
    #     decimals = max(0, oom - 6 + digits - 1)
    #     return @sprintf("")
    # end


    if x >= 100e6
        @sprintf("%.0fM", floor(x / 1e6))
    elseif x >= 100e3
        @sprintf("%.0fk", floor(x / 1e3))
    elseif x >= 1000
        @sprintf("%.0f", x)
    elseif x >= 100
        @sprintf("%.1f", x)
    elseif x >= 10
        @sprintf("%.2f", x)
    elseif x >= 1e-4
        @sprintf("%.4f", x)
    else
        @sprintf("%.2e", x)
    end
end

norm_sq(t) = dot(t, t)

# TODO move these to an Intervals module, together with all the other intervals stuff.

intersect_intervals(a, b) = (max(a[1], b[1]), min(a[2], b[2]))

is_subinterval(a, b) = a[1] >= b[1] && a[2] <= b[2]

end
