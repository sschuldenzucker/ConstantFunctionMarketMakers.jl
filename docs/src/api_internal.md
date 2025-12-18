# Internal API

```@meta
CurrentModule = ConstantFunctionMarketMakers
```

These functions are public but not exported. You can use them to define new AMM types or for specialized applications, but in general, the AMM / trading curve functions will be more convenient

```@docs
mk_tc
lt_vp
residual
```

The following functions query information from AMMs. They are supported but you usually want to create a [`TradingCurve`](@ref) instead and query that.

```@docs
t_l_p
l_t
p_l_t
```

## Other interface functions

```@docs
display_params
```

## AMM wrappers

```@docs
AMMWithMeta
TradingCurveWithMeta
ConcentratedAMM
```

