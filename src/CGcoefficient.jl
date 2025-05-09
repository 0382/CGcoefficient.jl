module CGcoefficient

export CG, CG0, threeJ, sixJ, nineJ, Racah,
       dCG, d3j, d6j, d9j, dRacah, xGaunt,
       lsjj, norm9J,
       Moshinsky,
       SqrtRational,
       exact_sqrt,
       simplify,
       iphase,
       is_same_parity,
       check_jm,
       check_couple,
       binomial_data_size, binomial_index,
       check_CG, check_3j, check_6j, check_9j,
       check_Gaunt, check_Moshinsky,
       fCG, fCG0, f3j, f6j, f9j, fRacah,
       fbinomial, unsafe_fbinomial,
       fMoshinsky, dfunc,
       flsjj, fCGspin, fCG3spin, fnorm9j,
       wigner_init_float


include("mpz.jl")
include("SqrtRational.jl")
include("util.jl")
include("WignerSymbols.jl")
include("floatWignerSymbols.jl")
include("PFRational/primetable.jl")
include("PFRational/PFRational.jl")
include("PFRational/pf_wigner.jl")
export PFRational, wigner_init_pf, pf_binomial, eCG, efCG, eCG0, efCG0, e3j, ef3j, e6j, ef6j

end # module CGcoefficient