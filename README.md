# ConstantFunctionMarketMakers.jl

Constant function market makers, as used on blockchains.

# Quickstart

This lets you plot some liquidity density curves for different AMMs:

```julia
using ConstantFunctionMarketMakers

amms = [
    ECLP(1-50e-4, 1+50e-4, 1-10e-4, 1000),
    C2LP(1-30e-4, 1+40e-4),
    StableSwap(400) |> concentrate(1-50e-4, 1+50e-4),
    CPMM(),
    CPSMM(0.9985),
  	CPSMM(0.9985) |> rate_scale(1.003, 1.0),
]
tcs = tc_v_p.(amms, 1e6, 1.0)

using Plots

p = plot()
for tc in tcs
  points = sample_dydp_p(tc; plim, unit=1e-4)
  plot!(p, points, label=short_str(tc; as_bp=true))
end
p
```

# Data Structures

There are two main abstract supertypes:

1. An `AMM` is the abstract configuration of a CFMM. Mathematically, this is just a function $F: \mathbb{R}^2 \to \mathbb{R}$ where the level sets are the trading curves, but computationally, it's not.

2. A `TradingCurve` is a specific trading curve of an AMM, together with an initial point on the curve. Mathematically, this is defined by an invariant $L$ plus a price $p$ and then the curve is equal to $F^(-1)(L)$, but computationally, there is more data.

All of these are immutable and composable. Specifically, you can take an AMM and run `concentrate(amm, alpha1, beta1)` to receive a new AMM that is concentrated into a sub-range of prices or `rate_scale(amm, rx, ry)` to receive an AMM that is rate-scaled by certain rates on its balances.

We support a number of calculations:

- Initialize trading curve of a given AMM based on balances or TVL and price.
- Query one balance given the other one (simulate trading), compute price at balance, balances at price, etc.
- Compute liquidity density and sample the liquidity density curve

Calculations are analytical where possible but some AMMs use numerical search for some operations.

# TODO / Not Implemented

- Support AMMs with more than two assets.
- Any kind of support for fees.
- Maybe support for mutable balances. Right now, everything is immutable and to adjust balances or rates, one creates a new trading curve. This could be inconvenient for larger iterations. (may not be worth it)

# License

MIT

# Development

We use the standard julia dependency management.

