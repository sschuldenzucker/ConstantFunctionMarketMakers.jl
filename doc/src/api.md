
# API Reference

```@meta
CurrentModule = ConstantFunctionMarketMakers
```

There are two main types of objects: [AMMs](#AMMs) and [Trading Curves](#Trading-Curves)

## AMMs

### AMM base type

The following is the supertype of all AMMs.

```@docs
AMM
```

### Base AMMs

```@docs
CPMM
C2LP
ECLP
StableSwap
WeightedPool
CPSMM
```

### Transformed AMMs

These take an AMM and yield an AMM that's transformed in some way.

```@docs
rate_scale
concentrate
named
flip
```

### Shorthands

```@docs
constant_rate_scaled_eclp
```

### Querying AMMs

AMMs themselves only store some basic information.

```@docs
alphabeta
kind_name
short_str
```

## Trading Curves

```@docs
TradingCurve
```

### Creating Trading Curves from AMMs

Given an AMM and parameters defining a starting point / invariant, one can create a `TradingCurve` from an AMM like so:

```@docs
tc_v_p
tc_l_p
tc_t
```

### Querying Trading Curves

With a trading curve constructed, we can query the curve in various ways.

```@docs
TCCommon
y_x
x_y
t_p
t_plus
v_init
```

## Helpers

```@docs
union_priceset
intersect_priceset
```

### Querying Liquidity Density

```@docs
dydp_p
dydlogp_p
dydp_max
dydlogp_max
sample_dydp_p
sample_dydlogp_p
```
