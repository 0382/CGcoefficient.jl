# Home

A package to calculate CG-coefficient and Wigner 3j, 6j, 9j symbols, and give exact results.

## Introduction

This package is inspired by Ref [^1]. See [CENS-MBPT](https://github.com/ManyBodyPhysics/CENS/blob/master/MBPT/VEffective/bhf-modules.f90) for details.

The idea is to simplify 3nj Symbols to sum combinations of binominal coefficients. We can calculate binominal coefficients by Pascal's Triangle, and store them first. Then we calculate 3nj Symbols using the stored binominal coefficients.

My code is just a toy model, using Julia's own `binominal` function, and use `BigInt` to get absolute results. If you want to get the exact result, and you don't use it for efficient numeric computing, it is really a good package.

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

In a markdown enviroment, such as jupyter notebook, it will give you a latex output. You can also do some calculations with the result.
```@example example
sixJ(1,2,3,4,5,6) * SqrtRational(1//7) * SqrtRational(1//13) * iphase(2+3+5+6)
```

In a console enviroment it will give out a text output.
```@repl example
nineJ(1,2,3,5,4,3,6,6,0)
```

You can also use `print` function to force print a text output.
```@example example
a = sixJ(1,2,3,2,1,2)
print(a)
```

## Index

- [Formula](formula.md)
- [API](api.md)

[^1]: T. Engeland and M. Hjorth-Jensen, the Oslo-FCI code. [https://github.com/ManyBodyPhysics/CENS](https://github.com/ManyBodyPhysics/CENS).