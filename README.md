# CGcoefficient.jl

A Julia package to calculate CG coefficient and Wigner 3j, 6j, 9j Symbols, and give exact results. To realize it, all the calculations are down by `BigInt` to avoid overflow, and it defines struct `SqrtRational` to store the exact results.


### Install

Just start a Julia REPL, and install it
```julia REPL
julia> ]
(@v1.5) pkg> add CGcoefficient
```

### Theory

[My blog](https://0382.github.io/2020/10/17/CG-coefficient-and-Wigner-3nj-Symbols/).

### Example

```julia REPL
julia> CG(1,2,3,1,1,2)
√2/√3
julia> nineJ(1,2,3,4,5,6,3,6,9)
√3/1274√5
```

More examples see: [jupyter example](./doc/example.ipynb)

### Reference

- [https://github.com/ManyBodyPhysics/CENS](https://github.com/ManyBodyPhysics/CENS)
- D. A. Varshalovich, A. N. Moskalev and V. K. Khersonskii, Quantum Theory of Angular Momentum, (World Scientific, 1988).