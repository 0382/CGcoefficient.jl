using Test
using CGcoefficient

include("test_sqrtrational.jl")
include("test_special.jl")
include("test_orthonormality.jl")
include("test_with_gsl.jl")
include("test_float_gsl.jl")
include("test_Moshinsky.jl")

@testset "test SqrtRational" begin test_sqrtrational() end
@testset "test SqrtRational show" begin test_show() end
@testset "test simplify" begin test_simplify() end

@testset "special condition: CG" begin test_special_CG(1//2:1//2:10) end
@testset "special condition: 6j" begin test_special_6j(1//2:1//2:5) end
@testset "special condition: Racah" begin test_special_Racah(1//2:1//2:5) end
@testset "special condition: 9j" begin test_special_9j(1//2:1//2:3) end

@testset "test orthonormality: dCG" begin test_orthonormality_dCG(0:8) end
@testset "test summation: d3j" begin test_summation_d3j(0:40) end
@testset "test orthonormality: d6j" begin test_orthonormality_d6j(0:6) end
@testset "test orthonormality: d9j" begin test_orthonormality_d9j(0:4) end

@testset "test Moshinsky" begin test_Moshinsky() end

try
    gsl3j(1, 1, 1, 0, 0, 0)
catch err
    if isa(err, ErrorException)
        println(err.msg)
        println("If you want to run test with libgsl, please download gsl library and rerun the test.")
        exit()
    end
    error(err)
end

@testset "test with gsl: 3j" begin test_3j_with_gsl(1:5) end
@testset "test with gsl: 6j" begin test_6j_with_gsl(1:5) end
@testset "test with gsl: 9j" begin test_9j_with_gsl(1:5) end
@testset "float version: f3j" begin test_f3j_with_gsl(50:55) end
@testset "float version: f6j" begin test_f6j_with_gsl(50:55) end
@testset "float version: f9j" begin test_f9j_with_gsl(50:52) end