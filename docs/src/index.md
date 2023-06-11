# Home

A package to calculate CG-coefficient, Racha coefficient, and Wigner 3j, 6j, 9j symbols. We offer three version API for all these coeddicients.

- **normal version**: `CG, threeJ, SixJ, Racah, nineJ`. Their parameters can be `Integer` or `Rational` (half integer), and their result is simplified by default.
- **double parameters version**: `dCG, d3j, d6j, dRacah d9j`. Their parameters means double of the real angular momentum, so can only be `Integer`. The result is not simplified. We will explain what is `simplify` later.
- **float version**: `fCG, f3j, f6j, fRacah, f9j`. Their parameters is same as double parameter version. They use `Float64` for calculation, so they are not exact but are fast for numeric calculation. Because these function use stored `binomial` result for speed up calculation, you should reserve space before calculate large angular momentum coefficients. See `wigner_init_float` function for details.

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
```@example example
d6j(2,4,6,8,10,12)
```

The `d` version functions do not simplify the result for the seek of speed, because `simplify` needs prime factorization which is slow. You can simplify the result explicitly,
```@example example
simplify(d6j(2,4,6,8,10,12))
```

You can also do some arithmetics with the result, thus do arithmetics using the `SqrtRational` type. The result is also not simplified
```@example example
x = sixJ(1,2,3,4,5,6) * exact_sqrt(1//7) * exact_sqrt(1//13) * iphase(2+3+5+6)
simplify(x)
```

In a console enviroment it will give out a text output.
```@repl example
nineJ(1,2,3,5,4,3,6,6,0)
```

You can also use `print` function to force print a text output.
```@example example
print(Racah(1,2,3,2,1,2))
```

## About

This package is inspired by Ref [^1]. See [CENS-MBPT](https://github.com/ManyBodyPhysics/CENS/blob/master/MBPT/VEffective/bhf-modules.f90) for details.

The idea is to simplify 3nj Symbols to sum combinations of binomial coefficients. We can calculate binomial coefficients by Pascal's Triangle, and store them first. Then we calculate 3nj Symbols using the stored binomial coefficients.

In this package, we just use the builtin `binomial` function for exact calculation. Only the float version uses stored `binomial`s.
## Index

- [Formula](formula.md)
- [API](api.md)

[^1]: T. Engeland and M. Hjorth-Jensen, the Oslo-FCI code. [https://github.com/ManyBodyPhysics/CENS](https://github.com/ManyBodyPhysics/CENS).