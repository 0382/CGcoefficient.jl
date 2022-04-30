# API

## Types
```@docs
HalfInt
SqrtRational
```

## Core functions
```@docs
CG
threeJ
sixJ
nineJ
Racah
```

People often use double of angular momentum quantum number as parameters, so we can use integer as parameters. This package also offers such functions, where the `d` letter means *double*.
```@docs
dCG
d3j
d6j
d9j
dRacah
```

## float version functions

Float version functions is always used for numeric calculation, so the parameters of all these functions (except `reserve_fbinomial`) are double of the exact angular momentum quantum number.
```@docs
fbinomial
fCG
f3j
f6j
f9j
fRacah
reserve_fbinomial
```

## Some useful function
```@docs
iphase
is_same_parity
check_jm
check_couple
exact_sqrt
simplify(::Integer)
simplify(::SqrtRational)
```