using Test
using BenchmarkTools

include("../src/CGcoefficient.jl")

const test_collect = 1:5

function gsl3j(dj1::Int, dj2::Int, dj3::Int, dm1::Int, dm2::Int, dm3::Int)
    ccall(
        (:gsl_sf_coupling_3j, "libgsl"),
        Cdouble,
        (Cint, Cint, Cint, Cint, Cint, Cint),
        dj1, dj2, dj3, dm1, dm2, dm3
    )
end

function gsl6j(dj1::Int, dj2::Int, dj3::Int, dj4::Int, dj5::Int, dj6::Int)
    ccall(
        (:gsl_sf_coupling_6j, "libgsl"),
        Cdouble,
        (Cint, Cint, Cint, Cint, Cint, Cint),
        dj1, dj2, dj3, dj4, dj5, dj6
    )
end

function gsl9j(dj1::Int, dj2::Int, dj3::Int,
                dj4::Int, dj5::Int, dj6::Int,
                dj7::Int, dj8::Int, dj9::Int)
    ccall(
        (:gsl_sf_coupling_9j, "libgsl"),
        Cdouble,
        (Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint),
        dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9
    )
end

function check_6j(dj1::Int, dj2::Int, dj3::Int, dj4::Int, dj5::Int, dj6::Int)
    check_couple(dj1, dj2, dj3) & check_couple(dj1, dj5, dj6) &
    check_couple(dj4, dj2, dj6) & check_couple(dj4, dj5, dj3)
end

function check_9j(dj1::Int, dj2::Int, dj3::Int, dj4::Int, dj5::Int, dj6::Int, dj7::Int, dj8::Int, dj9::Int)
    check_couple(dj1, dj2, dj3) & check_couple(dj4, dj5, dj6) & check_couple(dj7, dj8, dj9) &
    check_couple(dj1, dj4, dj7) & check_couple(dj2, dj5, dj8) & check_couple(dj3, dj6, dj9)
end

# 由于 gsl 库是浮点数运算，可能有误差，而我是准确计算。
# 有时我计算得到零，而 gsl 得到一个很小的数值，这种情况下 `≈` 运算会给出 false
# 这种情况应该排除
@testset "3j" begin
    for dj1 in test_collect
    for dj2 in test_collect
    for dj3 in test_collect
    for dm1 in -dj1:2:dj1
    for dm2 in -dj2:2:dj2
        dm3 = -dm1-dm2
        if check_couple(dj1, dj2, dj3) & check_jm(dj3, dm3)
            gsl = gsl3j(dj1, dj2, dj3, dm1, dm2, dm3)
            my = float(d3j(dj1, dj2, dj3, dm1, dm2, dm3))
            @test (my == 0) | (gsl ≈ my)
        end
    end end end end end
end

@testset "6j" begin
    for j1 in test_collect
    for j2 in test_collect
    for j3 in test_collect
    for j4 in test_collect
    for j5 in test_collect
    for j6 in test_collect
        if check_6j(j1, j2, j3, j4, j5, j6)
            gsl = gsl6j(j1, j2, j3, j4, j5, j6)
            my = float(d6j(j1,j2,j3,j4,j5,j6))
            @test (my == 0) | (gsl ≈ my)
        end
    end end end end end end
end

@testset "9j" begin
    for j1 in test_collect
    for j2 in test_collect
    for j3 in test_collect
    for j4 in test_collect
    for j5 in test_collect
    for j6 in test_collect
    for j7 in test_collect
    for j8 in test_collect
    for j9 in test_collect
        if check_9j(j1,j2,j3,j4,j5,j6,j7,j8,j9)
            gsl = gsl9j(j1,j2,j3,j4,j5,j6,j7,j8,j9)
            my = float(d9j(j1,j2,j3,j4,j5,j6,j7,j8,j9))
            @test (my == 0) | (gsl ≈ my)
        end
    end end end end end end end end end
end
