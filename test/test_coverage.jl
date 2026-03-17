# This file contains tests targeting code paths not covered by existing tests.

# Test utility functions directly
function test_utility_functions()
    # iphase
    @test iphase(0) == 1
    @test iphase(1) == -1
    @test iphase(2) == 1
    @test iphase(3) == -1

    # is_same_parity
    @test is_same_parity(0, 2) == true
    @test is_same_parity(1, 3) == true
    @test is_same_parity(0, 1) == false
    @test is_same_parity(2, 5) == false

    # check_jm
    @test check_jm(2, 2) == true
    @test check_jm(2, 0) == true
    @test check_jm(2, -2) == true
    @test check_jm(2, 4) == false   # |dm| > dj
    @test check_jm(2, 1) == false   # different parity

    # check_couple
    @test check_couple(2, 2, 2) == true
    @test check_couple(2, 2, 0) == true
    @test check_couple(2, 2, 6) == false  # triangle inequality violated
    @test check_couple(1, 1, 1) == false  # odd sum
    @test check_couple(0, 0, 2) == false  # triangle inequality

    # check_CG
    @test check_CG(2, 2, 2, 2, 0, 2) == true
    @test check_CG(2, 2, 2, 2, 0, 0) == false  # dm1+dm2 != dm3
    @test check_CG(2, 2, 6, 2, 2, 4) == false  # invalid couple

    # check_3j
    @test check_3j(2, 2, 2, 2, 0, -2) == true
    @test check_3j(2, 2, 2, 2, 0, 0) == false  # dm1+dm2+dm3 != 0
    @test check_3j(2, 2, 6, 2, 2, -4) == false  # invalid couple

    # check_6j
    @test check_6j(2, 2, 2, 2, 2, 2) == true
    @test check_6j(2, 2, 2, 2, 2, 10) == false  # invalid couple

    # check_9j
    @test check_9j(2, 2, 2, 2, 2, 2, 2, 2, 2) == true
    @test check_9j(2, 2, 2, 2, 2, 2, 2, 2, 10) == false

    # check_Gaunt
    @test check_Gaunt(1, 1, 2, 1, 0, -1) == true
    @test check_Gaunt(1, 1, 2, 2, 0, -2) == false  # |m1| > l1
    @test check_Gaunt(1, 1, 2, 1, 0, 0) == false   # m1+m2+m3 != 0
    @test check_Gaunt(1, 1, 5, 1, 0, -1) == false  # invalid couple

    # check_Moshinsky
    @test check_Moshinsky(0, 0, 0, 0, 0, 0, 0, 0, 0) == true
    @test check_Moshinsky(-1, 0, 0, 0, 0, 0, 0, 0, 0) == false  # N < 0
    @test check_Moshinsky(0, 0, -1, 0, 0, 0, 0, 0, 0) == false  # n < 0
    @test check_Moshinsky(0, 0, 0, 0, -1, 0, 0, 0, 0) == false  # n1 < 0
    @test check_Moshinsky(0, 0, 0, 0, 0, 0, -1, 0, 0) == false  # n2 < 0
    @test check_Moshinsky(0, 3, 0, 0, 0, 0, 0, 0, 0) == false   # Λ > L+l (triangle)
    @test check_Moshinsky(0, 0, 0, 0, 1, 0, 1, 0, 2) == false   # Λ > l1+l2 (triangle)
    @test check_Moshinsky(1, 0, 0, 0, 0, 0, 0, 0, 0) == false   # E+e != e1+e2
end

# Test binomial_data_size and binomial_index directly
function test_binomial_utils()
    # binomial_data_size
    @test binomial_data_size(0) == 1
    @test binomial_data_size(1) == 2
    @test binomial_data_size(2) == 4
    @test binomial_data_size(3) == 6
    @test binomial_data_size(4) == 9
    @test binomial_data_size(5) == 12
    @test binomial_data_size(6) == 16

    # binomial_index - the index for C(n, k) where k = min(k, n-k)
    @test binomial_index(0, 0) == 1
    @test binomial_index(1, 0) == 2
    # Let me just verify these are consistent with fbinomial
    wigner_init_float(10, "nmax", 0)
    for n in 0:10
        for k in 0:n
            idx = binomial_index(n, min(k, n-k))
            @test idx >= 1
            @test idx <= binomial_data_size(n)
        end
    end
end

# Test fbinomial edge cases
function test_fbinomial_edges()
    wigner_init_float(10, "nmax", 0)
    # valid cases
    @test fbinomial(10, 5) ≈ 252.0
    @test fbinomial(0, 0) ≈ 1.0
    @test fbinomial(5, 0) ≈ 1.0
    @test fbinomial(5, 5) ≈ 1.0
    # out-of-bounds cases should return 0
    @test fbinomial(-1, 0) == 0.0   # n < 0
    @test fbinomial(5, -1) == 0.0   # k < 0
    @test fbinomial(5, 6) == 0.0    # k > n
    @test fbinomial(10000, 5) == 0.0  # n >> any initialized nmax
end

# Test wigner_init_float all modes and error cases
function test_wigner_init_float_modes()
    # "Jmax" mode with different ranks
    wigner_init_float(5, "Jmax", 3)
    @test fCG(2, 2, 4, 2, 2, 4) ≈ 1.0
    wigner_init_float(5, "Jmax", 6)
    @test f6j(2, 2, 4, 2, 2, 2) ≈ float(d6j(2, 2, 4, 2, 2, 2))
    wigner_init_float(5, "Jmax", 9)
    @test f9j(2, 2, 2, 2, 2, 2, 2, 2, 2) ≈ float(d9j(2, 2, 2, 2, 2, 2, 2, 2, 2))
    # invalid rank for "Jmax" mode
    @test_throws ArgumentError wigner_init_float(5, "Jmax", 5)

    # "2bjmax" mode with different ranks
    wigner_init_float(10, "2bjmax", 3)
    @test fCG(2, 2, 4, 2, 2, 4) ≈ 1.0
    wigner_init_float(10, "2bjmax", 6)
    @test f6j(2, 2, 4, 2, 2, 2) ≈ float(d6j(2, 2, 4, 2, 2, 2))
    wigner_init_float(10, "2bjmax", 9)
    @test f9j(2, 2, 2, 2, 2, 2, 2, 2, 2) ≈ float(d9j(2, 2, 2, 2, 2, 2, 2, 2, 2))
    # invalid rank for "2bjmax" mode
    @test_throws ArgumentError wigner_init_float(5, "2bjmax", 5)

    # "nmax" mode
    wigner_init_float(20, "nmax", 0)
    @test fbinomial(20, 10) ≈ 184756.0

    # "Moshinsky" mode
    wigner_init_float(3, "Moshinsky", 0)
    @test fbinomial(19, 9) ≈ 92378.0

    # invalid mode
    @test_throws ArgumentError wigner_init_float(5, "invalid_mode", 3)
end

# Test wigner_init_pf all modes and error cases
function test_wigner_init_pf_modes()
    # "Jmax" mode with different ranks
    wigner_init_pf(3, "Jmax", 3)
    wigner_init_pf(3, "Jmax", 6)
    wigner_init_pf(3, "Jmax", 9)
    # invalid rank
    @test_throws ArgumentError wigner_init_pf(3, "Jmax", 5)

    # "2bjmax" mode
    wigner_init_pf(5, "2bjmax", 3)
    wigner_init_pf(5, "2bjmax", 6)
    wigner_init_pf(5, "2bjmax", 9)
    # invalid rank for "2bjmax" mode
    @test_throws ArgumentError wigner_init_pf(5, "2bjmax", 5)

    # "nmax" mode
    wigner_init_pf(10, "nmax", 0)
    @test numerator(pf_binomial(10, 5)) == 252

    # invalid mode
    @test_throws ArgumentError wigner_init_pf(5, "invalid_mode", 3)

    # pf_binomial error case: k < 0 always throws regardless of table size
    @test_throws ArgumentError pf_binomial(5, -1) # k < 0
end

# Test fCGspin with invalid inputs
function test_fCGspin_invalid()
    # invalid S (|S| > 1)
    @test fCGspin(1, 1, 2) == 0.0
    @test fCGspin(1, 1, -1) == 0.0  # unsigned(-1) > 1
    # invalid ds1
    @test fCGspin(2, 1, 1) == 0.0
    @test fCGspin(0, 1, 1) == 0.0
    # invalid ds2
    @test fCGspin(1, 2, 1) == 0.0
    @test fCGspin(1, 0, 1) == 0.0
end

# Test fCG3spin with invalid inputs
function test_fCG3spin_invalid()
    # S12 > 1
    @test fCG3spin(1, 1, 1, 2, 3) == 0.0
    # S12 == 0 but dS != 1
    @test fCG3spin(1, 1, 1, 0, 3) == 0.0
    # S12 == 1 but dS != 1 and dS != 3
    @test fCG3spin(1, 1, 1, 1, 5) == 0.0
    # invalid ds1, ds2, ds3
    @test fCG3spin(2, 1, 1, 1, 1) == 0.0
    @test fCG3spin(1, 2, 1, 1, 1) == 0.0
    @test fCG3spin(1, 1, 2, 1, 1) == 0.0
    @test fCG3spin(0, 1, 1, 1, 1) == 0.0
end

# Test SqrtRational error case: negative r
function test_sqrtrational_errors()
    # negative r should throw
    @test_throws ArgumentError SqrtRational(1//1, -1//1)
    # + operator with non-perfect-square r should throw
    x = SqrtRational(1//1, 3//1)  # √3
    @test_throws ArgumentError (x + 1)
end

# Test SqrtRational widen
function test_sqrtrational_widen()
    x = SqrtRational(3//5, 2//3)
    w = widen(x)
    @test float(x) ≈ float(w)
    @test w isa SqrtRational
end

# Test dfunc (Wigner d-function)
function test_dfunc()
    wigner_init_float(5, "Jmax", 9)

    # β = 0: d^j_{m1,m2}(0) = δ_{m1,m2}
    # For j=1/2 (dj=1): dfunc(1, dm1, dm2, 0.0) = δ_{dm1,dm2}
    @test dfunc(1, 1, 1, 0.0) ≈ 1.0
    @test dfunc(1, 1, -1, 0.0) ≈ 0.0
    @test dfunc(1, -1, 1, 0.0) ≈ 0.0
    @test dfunc(1, -1, -1, 0.0) ≈ 1.0

    # For j=1 (dj=2): dfunc(2, dm1, dm2, 0.0) = δ_{dm1,dm2}
    @test dfunc(2, 2, 2, 0.0) ≈ 1.0
    @test dfunc(2, 2, 0, 0.0) ≈ 0.0
    @test dfunc(2, 2, -2, 0.0) ≈ 0.0
    @test dfunc(2, 0, 0, 0.0) ≈ 1.0
    @test dfunc(2, -2, -2, 0.0) ≈ 1.0

    # For dj=0 (j=0): trivially 1
    @test dfunc(0, 0, 0, 0.0) ≈ 1.0
    arbitrary_angle = 1.234
    @test dfunc(0, 0, 0, arbitrary_angle) ≈ 1.0

    # Specific values for j=1/2 (dj=1):
    # dfunc uses β directly: dfunc(1, 1, 1, β) = cos(β)
    β = π/4
    @test dfunc(1, 1, 1, β) ≈ cos(β)
    @test dfunc(1, 1, -1, β) ≈ -sin(β)
    @test dfunc(1, -1, 1, β) ≈ sin(β)
    @test dfunc(1, -1, -1, β) ≈ cos(β)

    # Specific values for j=1 (dj=2):
    # dfunc(2, 2, 2, β) = cos²(β)
    β = π/6
    @test dfunc(2, 2, 2, β) ≈ cos(β)^2
    @test dfunc(2, 2, 0, β) ≈ -sin(2β)/sqrt(2)
    @test dfunc(2, 0, 0, β) ≈ cos(2β)

    # Invalid input (check_jm fails): returns 0
    @test dfunc(2, 3, 1, 0.5) == 0.0   # dm1 > dj
    @test dfunc(2, 1, 1, 0.5) == 0.0   # dm1 has wrong parity for dj=2

    # Symmetry: dfunc(dj, dm1, dm2, β) = dfunc(dj, -dm2, -dm1, β)
    β = 0.7
    for dj in [1, 2, 3, 4], dm1 in -dj:2:dj, dm2 in -dj:2:dj
        @test dfunc(dj, dm1, dm2, β) ≈ dfunc(dj, -dm2, -dm1, β)
    end
end

# Test Moshinsky with D != 1 (different mass ratio)
function test_moshinsky_D()
    wigner_init_float(5, "Moshinsky", 0)

    # Moshinsky bracket with D = 1 (same as default)
    m1 = Moshinsky(0, 0, 0, 0, 0, 0, 0, 0, 0, 1)
    @test m1 == 1

    # Moshinsky bracket with D = 2//1 (different mass ratio)
    # This exercises the D != 1 code path in _Moshinsky
    m2 = Moshinsky(0, 0, 0, 0, 0, 0, 0, 0, 0, 2//1)
    @test m2 == 1  # For n=l=N=L=n1=l1=n2=l2=0, Λ=0, the result is 1 regardless of D

    # A case with nonzero quantum numbers and D != 1
    m3 = Moshinsky(0, 1, 0, 0, 0, 0, 0, 1, 1, 2//1)
    @test m3 isa SqrtRational

    # Test that the float version also works
    fm3 = fMoshinsky(0, 1, 0, 0, 0, 0, 0, 1, 1, sqrt(2.0))
    @test fm3 isa Float64

    # Verify check: invalid Moshinsky returns 0 even with D parameter
    @test Moshinsky(-1, 0, 0, 0, 0, 0, 0, 0, 0, 2) == 0
end

# Test fRacah directly
function test_fRacah()
    wigner_init_float(5, "Jmax", 6)
    # fRacah is related to f6j by phase factor
    # Racah(j1, j2, j3, j4, j5, j6) = (-1)^{j1+j2+j3+j4} * sixJ(j1, j2, j5, j4, j3, j6)
    for j1 in (1, 2), j2 in (1, 2), j3 in (0, 1, 2), j4 in (1, 2), j5 in (1, 2), j6 in (0, 1, 2)
        if check_6j(j1, j2, j5, j4, j3, j6)
            exact = dRacah(j1, j2, j3, j4, j5, j6)
            fast = fRacah(j1, j2, j3, j4, j5, j6)
            @test fast ≈ float(exact)
        end
    end
end

# Test norm9J and fnorm9j directly
function test_norm9J_direct()
    # norm9J(j1,j2,j3, j4,j5,j6, j7,j8,j9) = sqrt((2j3+1)(2j6+1)(2j7+1)(2j8+1)) * nineJ(...)
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
