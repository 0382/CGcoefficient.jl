# Home

A package to calculate CG-coefficient, Racha coefficient, and Wigner 3j, 6j, 9j symbols. We offer two versions of API for all these coefficients.

1. The exact functions return `SqrtRational`, which are designed for demonstration. They use `BigInt` in the internal calculation, and do not cache the binomial table, so they are not efficient.
2. The floating-point functions return `Float64`, which are designed for numeric calculation. They use `Int, Float64` in the internal calculation, and you should pre-call [`wigner_init_float`](@ref) to calculate and cache the binomial table for later calculation. They may give inaccurate result for vary large angular momentum, due to floating-point arithmetic. You can find the max error at here: [wigner-benchmark](https://github.com/0382/wigner-benchmark).

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
In a markdown enviroment, such as jupyter notebook, it will give you a latex output.
You can also do some arithmetics with the result, thus do arithmetics using the `SqrtRational` type.
```@example example
x = sixJ(1,2,3,4,5,6) * exact_sqrt(1//7) * exact_sqrt(1//13) * iphase(2+3+5+6)
```

The result is not simplified by default, you can use `simplify` function to simplify it.
```@example example
simplify(x)
```

In a console enviroment it will give out a text output.
```@repl example
nineJ(1,2,3,5,4,3,6,6,0)
```

You can convert a `SqrtRational` in to `BigFloat`,
```@example example
t = sixJ(17/2,9,19/2,15/2,10,15/2)
float(t)
```

Calling `wigner_init_float` first to pre-calculate and cache binomial table.
```@example example
wigner_init_float(10, "Jmax", 6)
```
Then call a float version function
```@example example
f6j(17,18,19,15,20,15)
```

## About

This package is inspired by Ref [^1]. See [CENS-MBPT](https://github.com/ManyBodyPhysics/CENS/blob/master/MBPT/VEffective/bhf-modules.f90) for details.

The idea is to simplify 3nj Symbols to sum combinations of binomial coefficients. We can calculate binomial coefficients by Pascal's Triangle, and store them first. Then we calculate 3nj Symbols using the stored binomial coefficients.

In this package, we just use the builtin `binomial` function for exact calculation. Only the float version uses stored `binomial`s.

[^1]: T. Engeland and M. Hjorth-Jensen, the Oslo-FCI code. [https://github.com/ManyBodyPhysics/CENS](https://github.com/ManyBodyPhysics/CENS).
