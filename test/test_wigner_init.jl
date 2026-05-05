# this file tests wigner_init_float and wigner_init_pf:
# all supported modes ("Jmax", "2bjmax", "nmax", "Moshinsky") and error cases

function test_wigner_init_float_modes()
    # "Jmax" mode with all supported ranks
    wigner_init_float(5, "Jmax", 3)
    @test fCG(2, 2, 4, 2, 2, 4) ≈ 1.0
    wigner_init_float(5, "Jmax", 6)
    @test f6j(2, 2, 4, 2, 2, 2) ≈ float(d6j(2, 2, 4, 2, 2, 2))
    wigner_init_float(5, "Jmax", 9)
    @test f9j(2, 2, 2, 2, 2, 2, 2, 2, 2) ≈ float(d9j(2, 2, 2, 2, 2, 2, 2, 2, 2))
    # invalid rank for "Jmax" mode
    @test_throws ArgumentError wigner_init_float(5, "Jmax", 5)

    # "2bjmax" mode with all supported ranks
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

function test_wigner_init_pf_modes()
    # "Jmax" mode with all supported ranks
    wigner_init_pf(3, "Jmax", 3)
    wigner_init_pf(3, "Jmax", 6)
    wigner_init_pf(3, "Jmax", 9)
    # invalid rank for "Jmax" mode
    @test_throws ArgumentError wigner_init_pf(3, "Jmax", 5)

    # "2bjmax" mode with all supported ranks
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
    @test_throws ArgumentError pf_binomial(5, -1)
end
