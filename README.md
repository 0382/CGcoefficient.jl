# CGcoefficient.jl

A Julia package to calculate CG coefficient and Wigner 3j, 6j, 9j Symbols in an absolute way.

All the calculations are down by `BigInt` to avoid overflow. And define struct `SqrtRational` to store results.

### Theory

[My blog](https://0382.github.io/2020/10/17/CG-coefficient-and-Wigner-3nj-Symbols/).

### Example

```julia
julia> CG(1,2,3,1,1,2)
√2/√3
julia> nineJ(1,2,3,4,5,6,3,6,9)
√3/1274√5
```

More examples see: [jupyter example](./doc/example.ipynb)

### Reference

- [https://github.com/ManyBodyPhysics/CENS](https://github.com/ManyBodyPhysics/CENS)
- D. A. Varshalovich, A. N. Moskalev and V. K. Khersonskii, Quantum Theory of Angular Momentum, (World Scientific, 1988).