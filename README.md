# CGcoefficient.jl

[[中文](README_zh.md)]

A Julia package to calculate CG-coefficient and Wigner 3j, 6j and 9j symbols exactly. All the calculations are down by `BigInt` to avoid overflow, and defines struct `SqrtRational` to store the exact results.

For more details and the calculation formula, please see the [Document](https://0382.github.io/CGcoefficient.jl-docs/).

### Install

Just start a Julia REPL, and install it
```julia-repl
julia> ]
pkg> add CGcoefficient
```

### Example

```julia-repl
julia> CG(1,2,3,1,1,2)
√2/√3
julia> nineJ(1,2,3,4,5,6,3,6,9)
√3/1274√5
```

More examples see: [Document](https://0382.github.io/CGcoefficient.jl-docs/).

### Reference

- [https://github.com/ManyBodyPhysics/CENS](https://github.com/ManyBodyPhysics/CENS)
- D. A. Varshalovich, A. N. Moskalev and V. K. Khersonskii, Quantum Theory of Angular Momentum, (World Scientific, 1988).