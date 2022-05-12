push!(LOAD_PATH, "../src/")

using CGcoefficient
using BenchmarkTools
using Test

# wapper of gsl 3nj functions
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

function calculate_3j(func::Function, test_range::AbstractArray)
    if func == f3j
        reserve_fbinomial(cld(maximum(test_range), 2), "Jmax", 3)
    end
    sum = 0.0
    for dj1 in test_range
    for dj2 in test_range
    for dj3 in test_range
    for dm1 in -dj1:2:dj1
    for dm2 in -dj2:2:dj2
        dm3 = -dm1-dm2
        if check_couple(dj1, dj2, dj3) & check_jm(dj3, dm3)
            sum += func(dj1, dj2, dj3, dm1, dm2, dm3)
        end
    end end end end end
    return sum
end

function calculate_6j(func::Function, test_range::AbstractArray)
    if func == f6j
        reserve_fbinomial(cld(maximum(test_range), 2), "Jmax", 6)
    end
    sum = 0.0
    for j1 in test_range
    for j2 in test_range
    for j3 in test_range
    for j4 in test_range
    for j5 in test_range
    for j6 in test_range
        if check_6j(j1, j2, j3, j4, j5, j6)
            sum += func(j1, j2, j3, j4, j5, j6)
        end
    end end end end end end
    return sum
end

function calculate_9j(func::Function, test_range::AbstractArray)
    if func == f9j
        reserve_fbinomial(cld(maximum(test_range), 2), "Jmax", 9)
    end
    sum = 0.0
    for j1 in test_range
    for j2 in test_range
    for j3 in test_range
    for j4 in test_range
    for j5 in test_range
    for j6 in test_range
    for j7 in test_range
    for j8 in test_range
    for j9 in test_range
        if check_9j(j1,j2,j3,j4,j5,j6,j7,j8,j9)
            sum += func(j1,j2,j3,j4,j5,j6,j7,j8,j9)
        end
    end end end end end end end end end
    return sum
end

test_range = 40:50

t1 = @belapsed calculate_3j(f3j, test_range)
t2 = @belapsed calculate_3j(gsl3j, test_range)

println("diff = ", calculate_3j(f3j, test_range) - calculate_3j(gsl3j, test_range))
println("f3j time = $(t1)s")
println("gsl3j time = $(t2)s")

t1 = @belapsed calculate_6j(f6j, test_range)
t2 = @belapsed calculate_6j(gsl6j, test_range)

println("diff = ", calculate_6j(f6j, test_range) - calculate_6j(gsl6j, test_range))
println("f6j time = $(t1)s")
println("gsl6j time = $(t2)s")

test_range = 30:35

t1 = @belapsed calculate_9j(f9j, test_range)
t2 = @belapsed calculate_9j(gsl9j, test_range)

println("diff = ", calculate_9j(f9j, test_range) - calculate_9j(gsl9j, test_range))
println("f9j time = $(t1)s")
println("gsl9j time = $(t2)s")