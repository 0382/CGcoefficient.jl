# this file tests the Wigner d-function: dfunc(dj, dm1, dm2, β)

function test_dfunc()
    wigner_init_float(5, "Jmax", 9)

    # β = 0: d^j_{m1,m2}(0) = δ_{m1,m2}
    # For j=1/2 (dj=1)
    @test dfunc(1, 1, 1, 0.0) ≈ 1.0
    @test dfunc(1, 1, -1, 0.0) ≈ 0.0
    @test dfunc(1, -1, 1, 0.0) ≈ 0.0
    @test dfunc(1, -1, -1, 0.0) ≈ 1.0

    # For j=1 (dj=2)
    @test dfunc(2, 2, 2, 0.0) ≈ 1.0
    @test dfunc(2, 2, 0, 0.0) ≈ 0.0
    @test dfunc(2, 2, -2, 0.0) ≈ 0.0
    @test dfunc(2, 0, 0, 0.0) ≈ 1.0
    @test dfunc(2, -2, -2, 0.0) ≈ 1.0

    # For j=0 (dj=0): trivially 1 for any angle
    @test dfunc(0, 0, 0, 0.0) ≈ 1.0
    arbitrary_angle = 1.234
    @test dfunc(0, 0, 0, arbitrary_angle) ≈ 1.0

    # Specific analytic values for j=1/2 (dj=1):
    # d^{1/2}_{1/2,1/2}(β) = cos(β), d^{1/2}_{1/2,-1/2}(β) = -sin(β)
    β = π/4
    @test dfunc(1, 1, 1, β) ≈ cos(β)
    @test dfunc(1, 1, -1, β) ≈ -sin(β)
    @test dfunc(1, -1, 1, β) ≈ sin(β)
    @test dfunc(1, -1, -1, β) ≈ cos(β)

    # Specific analytic values for j=1 (dj=2):
    # d^1_{1,1}(β) = cos²(β), d^1_{0,0}(β) = cos(2β)
    β = π/6
    @test dfunc(2, 2, 2, β) ≈ cos(β)^2
    @test dfunc(2, 2, 0, β) ≈ -sin(2β)/sqrt(2)
    @test dfunc(2, 0, 0, β) ≈ cos(2β)

    # Invalid input (check_jm fails): returns 0.0
    @test dfunc(2, 3, 1, 0.5) == 0.0   # |dm1| > dj
    @test dfunc(2, 1, 1, 0.5) == 0.0   # dm1 has wrong parity for dj=2

    # Symmetry: d^j_{m1,m2}(β) = d^j_{-m2,-m1}(β)
    β = 0.7
    for dj in [1, 2, 3, 4], dm1 in -dj:2:dj, dm2 in -dj:2:dj
        @test dfunc(dj, dm1, dm2, β) ≈ dfunc(dj, -dm2, -dm1, β)
    end
end
