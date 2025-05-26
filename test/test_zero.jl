function test_trival_zero()
    x = 0.5
    @test CG(x, x, 1, x, x, -1) == 0
    @test CG(x, x, x, x, x, x) == 0
    @test CG(1, 2, 10, 0, 0, 0) == 0
    @test CG(-1, 2, 3, 0, 0, 0) == 0
    @test CG(x, x, 1, 1.5, x, 2) == 0
    @test CG(1, 2, 3, x, x, 1) == 0
    @test CG(2, 3, 4, 0, 0, 0) == 0

    @test sixJ(x, x, x, x, x, x) == 0
    @test sixJ(x, x, 1, x, x, x) == 0
    @test sixJ(x, x, 1, -1, x, 1) == 0
    @test sixJ(1, 1, 2, x, x, x) == 0

    @test nineJ(x, x, x, x, x, x, x, x, x) == 0
    @test nineJ(x, x, 1, x, x, x, x, x, x) == 0
    @test nineJ(x, x, 1, x, x, 1, x, x, x) == 0
    @test nineJ(x, x, 1, x, x, 1, x, x, 1) == 0
    @test nineJ(x, x, 1, x, x, 1, 1, x, x) == 0
    @test nineJ(x, x, 1, x, x, 0, 1, 1, 2) == 0

    @test xGaunt(1, 2, 3, 3, 2, 1) == 0
    @test xGaunt(1, 1, 1, 0, 0, 0) == 0

    @test Moshinsky(-1, 0, 0, 0, 0, 0, 0, 0, 0) == 0
    @test Moshinsky(3, 1, 8, 0, 5, 0, 6, 1, 9) == 0
    @test Moshinsky(3, 1, 4, 8, 5, 0, 6, 1, 9) == 0
    @test Moshinsky(3, 1, 4, 0, 5, 0, 6, 1, 1) == 0
end

function test_zero()
    @test sixJ(2, 2, 2, 1.5, 1.5, 1.5) == 0
    @test Racah(2, 2, 2, 2, 1, 3) == 0
    @test nineJ(3, 3, 3, 3, 3, 3, 3, 3, 3) == 0
end