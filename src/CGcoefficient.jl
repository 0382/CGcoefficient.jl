module CGcoefficient

export CG, CG0, threeJ, sixJ, nineJ, Racah,
       dCG, d3j, d6j, d9j, dRacah,
       lsjj, norm9J,
       Moshinsky,
       SqrtRational,
       exact_sqrt,
       simplify,
       HalfInt,
       iphase,
       is_same_parity,
       check_jm,
       check_couple,
       binomial_data_size, binomial_index,
       fCG, fCG0, f3j, f6j, f9j, fRacah,
       fbinomial, unsafe_fbinomial,
       fMoshinsky, dfunc,
       flsjj, fCGspin, fnorm9j,
       wigner_init_float

include("SqrtRational.jl")
include("util.jl")
include("WignerSymbols.jl")
include("floatWignerSymbols.jl")


let
    # Precompute at precompilation time
    _fbinomial_nmax[] = 67
    resize!(_fbinomial_data, binomial_data_size(get_fbinomial_nmax()))
    for n = 0:get_fbinomial_nmax()
        for k = 0:div(n, 2)
            get_fbinomial_data()[binomial_index(n, k)] = binomial(UInt64(n), UInt64(k))
        end
    end
end

end # module CGcoefficient
