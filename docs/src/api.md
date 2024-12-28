# API

## Types
```@docs
SqrtRational
```

## Some useful function
```@docs
iphase
is_same_parity
check_jm
check_couple
binomial_data_size
binomial_index
check_CG
check_3j
check_6j
check_9j
exact_sqrt
float(::SqrtRational)
simplify(::Integer)
simplify(::SqrtRational)
```

## Exact functions

The default functions give out exact result in the format of `SqrtRational`.
The results are simplified to give out shotest possible result.
Their arguments are `Real` (but only allow twice of which can be convert into integer).
```@docs
CG
CG0
threeJ
sixJ
nineJ
norm9J
lsjj
Racah
Moshinsky
```

People often use double of angular momentum quantum number as parameters, so we can use integer as parameters. This package also offers such functions, where the `d` letter means *double*.
These functions also give out exact `SqrtRational` results, but are not simplified.
Because the `simplify` function is quite slow, if you want to do some calculation for the result,
we suggest to use `d`-precedent functions first and `simplify` after call calculations.
```@docs
dCG
d3j
d6j
d9j
dRacah
```

## float version functions

Float version functions is always used for numeric calculation. They are designed for fast calculation.
You should call `wigner_init_float` to reserve the inner **binomial table**.
They only resive `Integer` arguments, thus `fCG, f3j, f6j, fRacha, f9j` only resive arguements which are
double of the exact quantum number. The rest functions do not need to do so.
The difference is labeled with the arguement name: `dj` means of double of the quantum number, while `j` means
the exact quantum number.
```@docs
fbinomial
unsafe_fbinomial
fCG
fCG0
fCGspin
fCG3spin
f3j
f6j
fRacah
f9j
fnorm9j
flsjj
fMoshinsky
dfunc
wigner_init_float
```

## prime factorization version

```@docs
PFRational
gcd(::PFRational, ::PFRational)
lcm(::PFRational, ::PFRational)
wigner_init_pf
pf_binomial
eCG
eCG0
e3j
e6j
efCG
efCG0
ef3j
ef6j
```
