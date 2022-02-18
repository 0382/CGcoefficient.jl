# this file test the behavior of `SqrtRational`

function test_show()
    io = IOBuffer()
    show(io, exact_sqrt(1))
    @test String(take!(io)) == "1"
    show(io, exact_sqrt(3//5))
    @test String(take!(io)) == "√(3//5)"
    show(io, "text/plain",  2//3 * exact_sqrt(5//6))
    @test String(take!(io)) == "(2//3)√(5//6)"
    show(io, "text/markdown", 2//3 * exact_sqrt(5//6))
    @test String(take!(io)) == "\$\\frac{2}{3}\\sqrt{\\frac{5}{6}}\$\n"
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