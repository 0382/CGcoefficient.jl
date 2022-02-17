# Home

A package to calculate CG-coefficient, Racha coefficient, and Wigner 3j, 6j, 9j symbols. It store the exact result with `SqrtRational` type.

## Introduction

This package is inspired by Ref [^1]. See [CENS-MBPT](https://github.com/ManyBodyPhysics/CENS/blob/master/MBPT/VEffective/bhf-modules.f90) for details.

The idea is to simplify 3nj Symbols to sum combinations of binominal coefficients. We can calculate binominal coefficients by Pascal's Triangle, and store them first. Then we calculate 3nj Symbols using the stored binominal coefficients. However, in current version, I just use the builtin `binominal` function, and calculation with large integer will overflow.

## Install

Just install with Julia REPL and enjoy it.

```julia-repl
pkg> add CGcoefficient
```

## Usage

```@setup example
push!(LOAD_PATH, "../../src/") # hide
```
```@example example
using CGcoefficient
sixJ(1,2,3,4,5,6)
```

In a markdown enviroment, such as jupyter notebook, it will give you a latex output. For the seek of speed, we do not simplify the result during the arithmetics. You can simplify the result explicitly.

```@example example
simplify(sixJ(1,2,3,4,5,6))
```

You can also do some arithmetics with the result, thus do arithmetics using the `SqrtRational` type.
```@example example
x = sixJ(1,2,3,4,5,6) * exact_sqrt(1//7) * exact_sqrt(1//13) * iphase(2+3+5+6)
simplify(x)
```

In a console enviroment it will give out a text output.
```@repl example
simplify(nineJ(1,2,3,5,4,3,6,6,0))
```

You can also use `print` function to force print a text output.
```@example example
a = Racah(1,2,3,2,1,2)
print(simplify(a))
```

## Index

- [Formula](formula.md)
- [API](api.md)

[^1]: T. Engeland and M. Hjorth-Jensen, the Oslo-FCI code. [https://github.com/ManyBodyPhysics/CENS](https://github.com/ManyBodyPhysics/CENS).