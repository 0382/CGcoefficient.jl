# this file test the validity by comparing with the widely used `libgsl`

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

# Because the gsl compute wigner 3nj symbols using float point number,
# sometimes my code give out zero, but gsl give out a very small number.
# In this condition, `≈` operator will give false result.
function test_3j_with_gsl(test_range::AbstractArray)
    for dj1 in test_range
    for dj2 in test_range
    for dj3 in test_range
    for dm1 in -dj1:2:dj1
    for dm2 in -dj2:2:dj2
        dm3 = -dm1-dm2
        if check_couple(dj1, dj2, dj3) & check_jm(dj3, dm3)
            gsl = gsl3j(dj1, dj2, dj3, dm1, dm2, dm3)
            my = float(d3j(dj1, dj2, dj3, dm1, dm2, dm3))
            @test (my == 0 && gsl < 1e-10) || (gsl ≈ my)
        end
    end end end end end
end

function test_6j_with_gsl(test_range::AbstractArray)
    for j1 in test_range
    for j2 in test_range
    for j3 in test_range
    for j4 in test_range
    for j5 in test_range
    for j6 in test_range
        if check_6j(j1, j2, j3, j4, j5, j6)
            gsl = gsl6j(j1, j2, j3, j4, j5, j6)
            my = float(d6j(j1,j2,j3,j4,j5,j6))
            @test (my == 0 && gsl < 1e-10) || (gsl ≈ my)
        end
    end end end end end end
end

function test_9j_with_gsl(test_range::AbstractArray)
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
            gsl = gsl9j(j1,j2,j3,j4,j5,j6,j7,j8,j9)
            my = float(d9j(j1,j2,j3,j4,j5,j6,j7,j8,j9))
            @test (my == 0 && gsl < 1e-10) | (gsl ≈ my)
        end
    end end end end end end end end end
end
