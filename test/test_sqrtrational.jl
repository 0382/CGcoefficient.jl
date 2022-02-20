# this file test the behavior of `SqrtRational`

function test_sqrtrational()
    x = 3//5 * exact_sqrt(2//3)
    @test zero(x) == 0
    @test 0*x == 0
    @test (0*x).r == 1
    @test (0*x) == zero(x)
    @test 1 == one(x)
    @test x*inv(x) == one(x)
    @test inv(x) == one(x)/x
    @test sign(+x) == 1
    @test sign(-x) == -1
    @test sign(0*x) == 0
    @test signbit(+x) == false
    @test signbit(-x) == true
    @test signbit(0*x) == false
end

function test_show()
    x = exact_sqrt(3)
    y = exact_sqrt(3//5)
    @test string(zero(x)) == "0"
    @test string(one(x)) == "1"
    @test string(2//3*one(x)) == "2//3"
    @test string(x) == "√3"
    @test string(y) == "√(3//5)"
    @test string(3*x) == "3√3"
    @test string(2//3*y) == "(2//3)√(3//5)"
    
    io = IOBuffer()
    show(io, "text/markdown", zero(x))
    @test String(take!(io)) == raw"$0$" * "\n"
    show(io, "text/markdown", one(x))
    @test String(take!(io)) == raw"$1$" * "\n"
    show(io, "text/markdown", 2//3*one(x))
    @test String(take!(io)) == raw"$\frac{2}{3}$" * "\n"
    show(io, "text/markdown", x)
    @test String(take!(io)) == raw"$\sqrt{3}$" * "\n"
    show(io, "text/markdown", y)
    @test String(take!(io)) == raw"$\sqrt{\frac{3}{5}}$" * "\n"
    show(io, "text/markdown", 3*x)
    @test String(take!(io)) == raw"$3\sqrt{3}$" * "\n"
    show(io, "text/markdown", 2//3*y)
    @test String(take!(io)) == raw"$\frac{2}{3}\sqrt{\frac{3}{5}}$" * "\n"
end

function test_simplify()
    @test simplify(2^7 * 3^3) == (2 * 3, 2^3 * 3)
    @test simplify(exact_sqrt(3//4)) == exact_sqrt(3) / 2
    test_simplify_with_6j()
end

function test_simplify_with_6j()
    test_range = 2:1//2:4
    for j1 in test_range
    for j2 in test_range
    for j3 in test_range
    for j4 in test_range
    for j5 in test_range
    for j6 in test_range
        try
            x = sixJ(j1,j2,j3,j4,j5,j6)
            @test simplify(x) == x
        catch err
            @test isa(err, ErrorException)
        end
    end end end end end end
end