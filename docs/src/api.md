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
check_Gaunt
check_Moshinsky
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
xGaunt
sixJ
nineJ
norm9J
lsjj
Racah
Moshinsky
```

People often use double of angular momentum quantum number as parameters, so we can use integer as parameters. This package also offers such functions, where the `d` letter means *double*.
These functions also give out exact `SqrtRational` results.
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
double of the exact quantum number. In some functions, e.g. `flsjj`, some quantum number can be half-integer,
the corresponding arguments start with `d` (which mean double of the quantum number);
while some other quantum number must be integer, the corresponding arguments does not starts with `d`.
```@docs
wigner_init_float
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
```

## prime factorization version

I implement a primary version of the prime factorization algorithm, based on the [wigxjpf](https://fy.chalmers.se/subatom/wigxjpf/) library.

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
