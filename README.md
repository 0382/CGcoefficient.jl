# CGcoefficient.jl

[[中文](README_zh.md)]

A package to calculate CG-coefficient, Racah coefficient, and Wigner 3j, 6j, 9j symbols. It store the exact result with `SqrtRational` type.

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
4√(1//24)

julia> simplify(ans)
√(2//3)

julia> simplify(nineJ(1,2,3,4,5,6,3,6,9))
1//1274√(3//5)
```

More examples see: [Document](https://0382.github.io/CGcoefficient.jl-docs/).

### Reference

- [https://github.com/ManyBodyPhysics/CENS](https://github.com/ManyBodyPhysics/CENS)
- D. A. Varshalovich, A. N. Moskalev and V. K. Khersonskii, Quantum Theory of Angular Momentum, (World Scientific, 1988).