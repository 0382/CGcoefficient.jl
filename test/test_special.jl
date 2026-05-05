# this file test the some special case if Wigner symbols
# especially the when one of arguement is equal to zero
# All formulat are refered in
# [1] D. A. Varshalovich, A. N. Moskalev and V. K. Khersonskii, Quantum Theory of Angular Momentum, (World Scientific, 1988).

# test special condition for CG coefficients
# Ref[1], P248, Sec 8.5, Formula(1)
function test_special_CG(test_range::AbstractArray)
    for j = test_range
        for m = -j:1:j
            cg = CG(j, j, 0, m, -m, 0)
            mycg = iphase(Int(j - m)) / exact_sqrt(2j + 1)
            @test cg == mycg
        end
    end
end

function test_CG0(test_range::AbstractArray)
    for j1 = test_range, j2 = test_range, j3 = test_range
        if check_couple(Int(2j1), Int(2j2), Int(2j3))
            cg0 = CG0(j1, j2, j3)
            tj = threeJ(j1, j2, j3, 0, 0, 0) * iphase(2j1 + j3) * exact_sqrt(2j3 + 1)
            @test cg0 == tj
        end
    end
end

function test_CGspin()
    wigner_init_float(2, "Jmax", 3)
    @test fCGspin(1, 1, 1) ≈ fCG(1, 1, 2, 1, 1, 2) ≈ 1.0
    @test fCGspin(1, 1, 0) ≈ fCG(1, 1, 0, 1, 1, 2) ≈ 0.0
    @test fCGspin(1, -1, 1) ≈ fCG(1, 1, 2, 1, -1, 0) ≈ 1/√2
    @test fCGspin(1, -1, 0) ≈ fCG(1, 1, 0, 1, -1, 0) ≈ 1/√2
    @test fCGspin(-1, 1, 1) ≈ fCG(1, 1, 2, -1, 1, 0) ≈ 1/√2
    @test fCGspin(-1, 1, 0) ≈ fCG(1, 1, 0, -1, 1, 0) ≈ -1/√2
    @test fCGspin(-1, -1, 1) ≈ fCG(1, 1, 2, -1, -1, -2) ≈ 1.0
    @test fCGspin(-1, -1, 0) ≈ fCG(1, 1, 0, -1, -1, -2) ≈ 0.0
end

function test_CG3spin()
    wigner_init_float(3, "Jmax", 3)
    for ds1 = (-1, 1) for ds2 = (-1, 1) for ds3 = (-1, 1)
        for S12 = (0,1)
            for dS = (1,3)
                x = fCG3spin(ds1, ds2, ds3, S12, dS)
                y = fCG(1, 1, 2S12, ds1, ds2, ds1+ds2) * fCG(2S12, 1, dS, ds1+ds2, ds3, ds1+ds2+ds3)
                @test x ≈ y
            end
        end
    end end end
end

function test_3j_as_CG(test_range::AbstractArray)
    for dj1 in test_range, dj2 in test_range, dj3 in test_range
        for dm1 in -dj1:2:dj1, dm2 in -dj2:2:dj2
            dm3 = -dm1 - dm2
            if check_couple(dj1, dj2, dj3) & check_jm(dj3, dm3)
                tj = d3j(dj1, dj2, dj3, dm1, dm2, dm3)
                cg = dCG(dj1, dj2, dj3, -dm1, -dm2, dm3)
                @test tj == cg * SqrtRational(iphase(dj1 + div(dj3 + dm3, 2)), 1//(dj3 + 1))
            end
        end
    end
    wigner_init_float(cld(maximum(test_range), 2), "Jmax", 3)
    for dj1 in test_range, dj2 in test_range, dj3 in test_range
        for dm1 in -dj1:2:dj1, dm2 in -dj2:2:dj2
            dm3 = -dm1 - dm2
            if check_couple(dj1, dj2, dj3) & check_jm(dj3, dm3)
                tj = f3j(dj1, dj2, dj3, dm1, dm2, dm3)
                cg = fCG(dj1, dj2, dj3, -dm1, -dm2, dm3)
                @test tj ≈ cg * iphase(dj1 + div(dj3 + dm3, 2)) / sqrt(dj3 + 1)
            end
        end
    end
end

# test special condition for 6j symbols
# Ref[1], P299, Sec 9.5, Formula (1)
function test_special_6j(test_range::AbstractArray)
    for j1 in test_range, j2 in test_range, j3 in test_range
        if check_couple(Int(2j1), Int(2j2), Int(2j3))
            sj = sixJ(j1, j2, j3, j2, j1, 0)
            mysj = iphase(Int(j1 + j2 + j3)) / exact_sqrt((2j1 + 1) * (2j2 + 1))
            @test sj == mysj
        end
    end
end

# test special condition for Racah coefficients
# Ref[1], P300, Sec 9.5, Formula (2)
function test_special_Racah(test_range::AbstractArray)
    for j1 in test_range, j2 in test_range, j3 in test_range
        if check_couple(Int(2j1), Int(2j2), Int(2j3))
            rc = Racah(0, j1, j2, j3, j1, j2)
            myrc = 1 / exact_sqrt((2j1 + 1) * (2j2 + 1))
            @test rc == myrc
        end
    end
end

# test special condition for 9j symbols
# Ref[1], P357, Sec 10.9, Formula (1)
function test_special_9j(test_range::AbstractArray)
    for j1 in test_range, j2 in test_range, j3 in test_range,
        j4 in test_range, j5 in test_range, j7 in test_range

        if check_6j(convert(Int, 2j1), convert(Int, 2j2), convert(Int, 2j3), convert(Int, 2j5), convert(Int, 2j4), convert(Int, 2j7))
            nj = nineJ(j1, j2, j3, j4, j5, j3, j7, j7, 0)
            snj = iphase(Int(j2 + j3 + j4 + j7)) / exact_sqrt((2j3 + 1) * (2j7 + 1))
            snj *= sixJ(j1, j2, j3, j5, j4, j7)
            @test nj == snj
        end
    end
end

function test_lsjj(test_range::AbstractArray)
    for l1 in test_range
        for l2 in test_range
            for L in abs(l1-l2):(l1+l2)
                for (j1, j2) in [(l1+1//2, l2+1//2), (l1+1//2, l2-1//2), (l1-1//2, l2+1//2), (l1-1//2, l2-1//2)]
                    for (S, J) in [(0, L), (1, L-1), (1, L), (1, L+1)]
                        @test lsjj(l1, l2, j1, j2, L, S, J) == norm9J(l1, 1//2, j1, l2, 1//2, j2, L, S, J)
                    end
                end
            end
        end
    end
end

function test_flsjj(test_range::AbstractArray)
    wigner_init_float(maximum(test_range), "Jmax", 9)
    for l1 in test_range
        for l2 in test_range
            for L in abs(l1-l2):(l1+l2)
                for (dj1, dj2) in [(2l1+1, 2l2+1), (2l1+1, 2l2-1), (2l1-1, 2l2+1), (2l1-1, 2l2-1)]
                    for (S, J) in [(0, L), (1, L-1), (1, L), (1, L+1)]
                        @test flsjj(l1, l2, dj1, dj2, L, S, J) ≈ fnorm9j(2l1, 1, dj1, 2l2, 1, dj2, 2L, 2S, 2J)
                    end
                end
            end
        end
    end
end

# Test fCGspin with invalid inputs (should return 0.0)
function test_fCGspin_invalid()
    # invalid S (|S| > 1)
    @test fCGspin(1, 1, 2) == 0.0
    @test fCGspin(1, 1, -1) == 0.0
    # invalid ds1 (must be ±1)
    @test fCGspin(2, 1, 1) == 0.0
    @test fCGspin(0, 1, 1) == 0.0
    # invalid ds2 (must be ±1)
    @test fCGspin(1, 2, 1) == 0.0
    @test fCGspin(1, 0, 1) == 0.0
end

# Test fCG3spin with invalid inputs (should return 0.0)
function test_fCG3spin_invalid()
    # S12 > 1
    @test fCG3spin(1, 1, 1, 2, 3) == 0.0
    # S12 == 0 but dS != 1
    @test fCG3spin(1, 1, 1, 0, 3) == 0.0
    # S12 == 1 but dS is not 1 or 3
    @test fCG3spin(1, 1, 1, 1, 5) == 0.0
    # invalid ds1, ds2, ds3 (must be ±1)
    @test fCG3spin(2, 1, 1, 1, 1) == 0.0
    @test fCG3spin(1, 2, 1, 1, 1) == 0.0
    @test fCG3spin(1, 1, 2, 1, 1) == 0.0
    @test fCG3spin(0, 1, 1, 1, 1) == 0.0
end

# Test fRacah directly against the exact dRacah
function test_fRacah()
    wigner_init_float(5, "Jmax", 6)
    for j1 in (1, 2), j2 in (1, 2), j3 in (0, 1, 2), j4 in (1, 2), j5 in (1, 2), j6 in (0, 1, 2)
        if check_6j(j1, j2, j5, j4, j3, j6)
            @test fRacah(j1, j2, j3, j4, j5, j6) ≈ float(dRacah(j1, j2, j3, j4, j5, j6))
        end
    end
end

# Test norm9J against its definition: sqrt((2j3+1)(2j6+1)(2j7+1)(2j8+1)) * nineJ(...)
function test_norm9J_direct()
    for j1 in (0, 1//2, 1), j2 in (0, 1//2, 1), j4 in (0, 1//2, 1), j5 in (0, 1//2, 1)
        for j3 in abs(j1-j2):(j1+j2), j6 in abs(j4-j5):(j4+j5)
            for j7 in abs(j1-j4):(j1+j4), j8 in abs(j2-j5):(j2+j5)
                for j9 in abs(j7-j8):(j7+j8)
                    if check_9j(Int(2j1),Int(2j2),Int(2j3),Int(2j4),Int(2j5),Int(2j6),Int(2j7),Int(2j8),Int(2j9)) &&
                       check_couple(Int(2j3),Int(2j6),Int(2j9))
                        n9 = norm9J(j1,j2,j3,j4,j5,j6,j7,j8,j9)
                        expected = exact_sqrt((2j3+1)*(2j6+1)*(2j7+1)*(2j8+1)) * nineJ(j1,j2,j3,j4,j5,j6,j7,j8,j9)
                        @test n9 == expected
                    end
                end
            end
        end
    end
end