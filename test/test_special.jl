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


function check_6j(j1::HalfInt, j2::HalfInt, j3::HalfInt,
    j4::HalfInt, j5::HalfInt, j6::HalfInt)
    dj1, dj2, dj3, dj4, dj5, dj6 = Int64.((2j1, 2j2, 2j3, 2j4, 2j5, 2j6))
    check_couple(dj1, dj2, dj3) & check_couple(dj1, dj5, dj6) &
    check_couple(dj4, dj2, dj6) & check_couple(dj4, dj5, dj3)
end

# test special condition for 9j symbols
# Ref[1], P357, Sec 10.9, Formula (1)
function test_special_9j(test_range::AbstractArray)
    for j1 in test_range, j2 in test_range, j3 in test_range,
        j4 in test_range, j5 in test_range, j7 in test_range

        if check_6j(j1, j2, j3, j5, j4, j7)
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