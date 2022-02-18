using Test
using CGcoefficient

include("test_sqrtrational.jl")
include("test_special.jl")
include("test_with_gsl.jl")

@testset "test SqrtRational show" begin test_show() end
@testset "test simplify" begin test_simplify() end

@testset "special condition: CG" begin test_special_CG(1//2:1//2:10) end
@testset "special condition: 3j" begin test_special_3j(1//2:1//2:10) end
@testset "special condition: 6j" begin test_special_6j(1//2:1//2:5) end
@testset "special condition: 9j" begin test_special_9j(1//2:1//2:5) end

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
@testset "test_with gsl: 6j" begin test_6j_with_gsl(1:5) end
@testset "test_with gsl: 9j" begin test_9j_with_gsl(1:5) end