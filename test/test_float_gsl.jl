# this file test the validity by comparing with the widely used `libgsl`

function test_f3j_with_gsl(test_range::AbstractArray)
    wigner_init_float(cld(maximum(test_range), 2), "Jmax", 3)
    for dj1 in test_range
    for dj2 in test_range
    for dj3 in test_range
    for dm1 in -dj1:2:dj1
    for dm2 in -dj2:2:dj2
        dm3 = -dm1-dm2
        if check_couple(dj1, dj2, dj3) & check_jm(dj3, dm3)
            gsl = gsl3j(dj1, dj2, dj3, dm1, dm2, dm3)
            my = f3j(dj1, dj2, dj3, dm1, dm2, dm3)
            @test isapprox(my, gsl; atol=1e-10)
        end
    end end end end end
end

function test_f6j_with_gsl(test_range::AbstractArray)
    wigner_init_float(cld(maximum(test_range), 2), "Jmax", 6)
    for j1 in test_range
    for j2 in test_range
    for j3 in test_range
    for j4 in test_range
    for j5 in test_range
    for j6 in test_range
        if check_6j(j1, j2, j3, j4, j5, j6)
            gsl = gsl6j(j1, j2, j3, j4, j5, j6)
            my = f6j(j1,j2,j3,j4,j5,j6)
            @test isapprox(my, gsl; atol=1e-10)
        end
    end end end end end end
end

function test_f9j_with_gsl(test_range::AbstractArray)
    wigner_init_float(cld(maximum(test_range), 2), "Jmax", 9)
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
            my = f9j(j1,j2,j3,j4,j5,j6,j7,j8,j9)
            @test isapprox(my, gsl; atol=1e-10)
        end
    end end end end end end end end end
end
