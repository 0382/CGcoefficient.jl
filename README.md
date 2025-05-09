# CGcoefficient.jl

[![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](LICENSE)
[![CI](https://github.com/0382/CGcoefficient.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/0382/CGcoefficient.jl/actions/workflows/CI.yml)
[![codecov.io](https://codecov.io/gh/0382/CGcoefficient.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/0382/CGcoefficient.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://0382.github.io/CGcoefficient.jl/dev)

[[中文](README_zh.md)]

A package to calculate CG-coefficient, Racah coefficient, Wigner 3j, 6j, 9j symbols and Moshinsky brakets.

One can get the exact result with `SqrtRational` type, which use `BigInt` to avoid overflow. And we also offer float version for numeric calculation, which is several times faster than [GNU Scientific Library](https://www.gnu.org/software/gsl/).

I also rewrite the float version with c++ for numeric calculation: [WignerSymbol](https://github.com/0382/WignerSymbol).

For more details and the calculation formula, please see the document [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://0382.github.io/CGcoefficient.jl/dev).

### Install

Just start a Julia REPL, and install it
```julia-repl
julia> ]
pkg> add CGcoefficient
```

### Example

```julia-repl
julia> CG(1,2,3,1,1,2)
√(2//3)

julia> nineJ(1,2,3,4,5,6,3,6,9)
1//1274√(3//5)

julia> f6j(6,6,6,6,6,6)
-0.07142857142857142
```

For more examples please see the document.

### API

This package contains two types of functions:

1. The exact functions return `SqrtRational`, which are designed for demonstration. They use `BigInt` in the internal calculation, and do not cache the binomial table, so they are not efficient.
2. The floating-point functions return `Float64`, which are designed for numeric calculation. They use `Int, Float64` in the internal calculation, and you should pre-call [`wigner_init_float`](https://0382.github.io/CGcoefficient.jl/stable/api/#CGcoefficient.wigner_init_float) to calculate and cache the binomial table for later calculation. They may give inaccurate result for vary large angular momentum, due to floating-point arithmetic. You can find the max error at here: [wigner-benchmark](https://github.com/0382/wigner-benchmark).

#### Exact functions

Exact function only for demonstration, some of them can be `Real` arguments, like `1, 3//2, 0.5`, but `4//3, 0.6` will throw error.

- `CG(j1, j2, j3, m1, m2, m3)`, CG-coefficient, `Real` arguments.
- `CG0(j1, j2, j3)`, CG-coefficient for `m1 = m2 = m3 = 0`, only integer angular momentum number is meaningful.
- `threeJ(j1, j2, j3, m1, m2, m3)`, Wigner 3j-symbol, `Real` arguments.
- `sixJ(j1, j2, j3, j4, j5, j6)`, Wigner 6j-symbol, `Real` arguments.
- `Racah(j1, j2, j3, j4, j5, j6)`, Racah coefficient, `Real` arguments.
- `nineJ(j1, j2, j3, j4, j5, j6, j7, j8, j9)`, Wigner 9j-symbol, `Real` arguments.
- `norm9J(j2, j3, j4, j5, j5, j6, j7, j8, j9)`, normalized 9j-symbol, `Real` arguments.
- `lsjj(l1, l2, j1, j2, L, S, J)`, LS-coupling to jj-coupling transform coefficient. It actually equals to a normalized 9j-symbol, but easy to use and faster. `j1, j2` can be `Real`.
- `Moshinsky(N, L, n, l, n1, l1, n2, l2, Λ, D)`, Moshinsky brakets, `Integer` arguments, except `D` is a rational (default is `1`).

#### Float functions

For faster numeric calculation, if the angular momentum number can be half-integer, the argument of the functions is actually double of the number. In this way, all arguments can be integers. In this package, the doubled arguments are named starts with `d`.

- `fCG(dj1, dj2, dj3, dm1, dm2, dm3)`, CG-coefficient.
- `fCG0(j1, j2, j3)`, CG-coefficient for `m1 = m2 = m3 = 0`.
- `fCGspin(ds1, ds2, S)`, quicker CG-coefficient for two spin-1/2 coupling.
- `fCG3spin(ds1, ds2, ds3, S12, dS)`, quickly evaluate `<S12,M12|1/2,m1;1/2,m2><S,M|S12,M12;1/2,m3>`.
- `f3j(dj1, dj2, dj3, dm1, dm2, dm3)`, Wigner 3j-symbol.
- `f6j(dj1, dj2, dj3, dj4, dj5, dj6)`, Wigner 6j-symbol.
- `fRacah(dj1, dj2, dj3, dj4, dj5, dj6)`, Racah coefficient.
- `f9j(dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9)`, Wigner 9j-symbol.
- `fnorm9j(dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9)`, normalized 9j-symbol.
- `flsjj(l1, l2, dj1, dj2, L, S, J)`, LS-coupling to jj-coupling transform coefficient.
- `fMoshinsky(N, L, n, l, n1, l1, n2, l2, Λ)`, Moshinsky brakets.
- `dfunc(dj, dm1, dm2, β)`, Wigner d-function.

### Reference

- [https://github.com/ManyBodyPhysics/CENS](https://github.com/ManyBodyPhysics/CENS)
- L. Wei, Unified approach for exact calculation of angular momentum coupling and recoupling coefficients, Computer Physics Communications 120, 222 (1999).
- L. Wei, New formula for 9-j symbols and their direct calculation, Computers in Physics 12, 632 (1998).
- D. A. Varshalovich, A. N. Moskalev and V. K. Khersonskii, *Quantum Theory of Angular Momentum*, (World Scientific, 1988).
- Buck et al. Nuc. Phys. A 600 (1996) 387-402.
- H. T. Johansson and C. Forssén, Fast and Accurate Evaluation of Wigner 3$j$, 6$j$, and 9$j$ Symbols Using Prime Factorization and Multiword Integer Arithmetic, SIAM J. Sci. Comput. 38, A376 (2016).
