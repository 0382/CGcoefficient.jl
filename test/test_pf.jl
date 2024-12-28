function test_PFRational()
    wigner_init_pf(10, "nmax", 0)
    x = pf_binomial(10, 5) # 252
    @test numerator(x) == 252
    @test denominator(x) == 1
    @test convert(Rational, x) == 252//1
    y = pf_binomial(10, 6) # 210
    @test numerator(gcd(x, y)) == 42
    @test numerator(lcm(x, y)) == 1260
    @test numerator(x * y) == 52920
    @test numerator(x / y) == 6
    @test denominator(x / y) == 5
    @test numerator(inv(x)) == 1
    @test denominator(inv(x)) == 252
    @test numerator(x^2) == 63504
end

function test_eCG(test_range::AbstractArray)
    wigner_init_pf(cld(maximum(test_range), 2), "Jmax", 3)
    for dj1 in test_range
    for dj2 in test_range
    for dj3 in test_range
    for dm1 in -dj1:2:dj1
    for dm2 in -dj2:2:dj2
        dm3 = dm1 + dm2
        d = dCG(dj1, dj2, dj3, dm1, dm2, dm3)
        e = eCG(dj1, dj2, dj3, dm1, dm2, dm3)
        ef = efCG(dj1, dj2, dj3, dm1, dm2, dm3)
        d = simplify(d)
        @test d == e
        @test convert(Float64, float(d)) == ef
    end end end end end
end

function test_eCG0(test_range::AbstractArray)
    wigner_init_pf(maximum(test_range), "Jmax", 3)
    for j1 in test_range
    for j2 in test_range
    for j3 in test_range
        d = CG0(j1, j2, j3)
        e = eCG0(j1, j2, j3)
        ef = efCG0(j1, j2, j3)
        @test d == e
        @test convert(Float64, float(d)) == ef
    end end end
end

function test_e3j(test_range::AbstractArray)
    wigner_init_pf(cld(maximum(test_range), 2), "Jmax", 3)
    for dj1 in test_range
    for dj2 in test_range
    for dj3 in test_range
    for dm1 in -dj1:2:dj1
    for dm2 in -dj2:2:dj2
        dm3 = -dm1-dm2
        d = d3j(dj1, dj2, dj3, dm1, dm2, dm3)
        e = e3j(dj1, dj2, dj3, dm1, dm2, dm3)
        ef = ef3j(dj1, dj2, dj3, dm1, dm2, dm3)
        d = simplify(d)
        @test d == e
        @test convert(Float64, float(d)) == ef
    end end end end end
end

function test_e6j(test_range::AbstractArray)
    wigner_init_pf(cld(maximum(test_range), 2), "Jmax", 6)
    for j1 in test_range
    for j2 in test_range
    for j3 in test_range
    for j4 in test_range
    for j5 in test_range
    for j6 in test_range
        d = d6j(j1, j2, j3, j4, j5, j6)
        e = e6j(j1, j2, j3, j4, j5, j6)
        ef = ef6j(j1, j2, j3, j4, j5, j6)
        d = simplify(d)
        if d != e
            println(j1, j2, j3, j4, j5, j6)
            println(d)
            println(e)
            return
        end
        @test d == e
        @test convert(Float64, float(d)) == ef
    end end end end end end
end


function test_big_binomial()
    wigner_init_pf(100, "nmax", 0)
    wigner_init_float(100, "nmax", 0)
    n = 100
    for k in 0:n
        x = pf_binomial(n, k)
        y = binomial(big(n), big(k))
        @test numerator(x) == y
        f = fbinomial(n, k)
        @test convert(Float64, numerator(x)) ≈ f
    end
end