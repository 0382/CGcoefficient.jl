# this file tests utility functions: iphase, parity checks, angular momentum checks,
# and binomial coefficient utilities (binomial_data_size, binomial_index, fbinomial)

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

# Test binomial_data_size and binomial_index
function test_binomial_utils()
    # binomial_data_size(n) counts stored entries up to C(n, *)
    @test binomial_data_size(0) == 1
    @test binomial_data_size(1) == 2
    @test binomial_data_size(2) == 4
    @test binomial_data_size(3) == 6
    @test binomial_data_size(4) == 9
    @test binomial_data_size(5) == 12
    @test binomial_data_size(6) == 16

    # binomial_index should stay within the bounds of binomial_data_size
    wigner_init_float(10, "nmax", 0)
    for n in 0:10
        for k in 0:n
            idx = binomial_index(n, min(k, n - k))
            @test idx >= 1
            @test idx <= binomial_data_size(n)
        end
    end
end

# Test fbinomial edge cases (out-of-range inputs should return 0.0)
function test_fbinomial_edges()
    wigner_init_float(10, "nmax", 0)
    # valid cases
    @test fbinomial(10, 5) ≈ 252.0
    @test fbinomial(0, 0) ≈ 1.0
    @test fbinomial(5, 0) ≈ 1.0
    @test fbinomial(5, 5) ≈ 1.0
    # out-of-bounds: should return 0.0
    @test fbinomial(-1, 0) == 0.0   # n < 0
    @test fbinomial(5, -1) == 0.0   # k < 0
    @test fbinomial(5, 6) == 0.0    # k > n
    @test fbinomial(10000, 5) == 0.0  # n far exceeds initialized nmax
end
