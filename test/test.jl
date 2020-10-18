using Test

include("../src/CGcoefficient.jl")

const test_collect1 = 1//2:1//2:10
const test_collect2 = 1//2:1//2:5


# 测试 CG 系数的特殊情况
@testset "CG" begin
for j = test_collect1
    for m = -j:1:j
        cg = CG(j, j, 0, m, -m, 0)
        mycg = iphase(Int(j-m)) / SqrtRational(2j+1)
        @test cg == mycg
    end
end
end

# 测试 3j 系数的特殊情况
@testset "3j" begin
for j = test_collect1
    for m = -j:1:j
        tj = threeJ(j, 1, j, -m, 0, m)
        mytj = iphase(Int(j-m)) * m / SqrtRational(j*(2j+1)*(j+1))
        @test tj == mytj
    end
end
end

# 测试 6j 系数的特殊情况
@testset "6j" begin
for j1 = test_collect2, j2 = test_collect2, j3 = test_collect2
    if check_couple(Int(2j1), Int(2j2), Int(2j3))
        sj = sixJ(j1, j2, j3, j2, j1, 0)
        mysj = iphase(Int(j1+j2+j3)) / SqrtRational((2j1+1)*(2j2+1))
        @test sj == mysj
    end
end
end


function check_6j(j1::AngularType, j2::AngularType, j3::AngularType,
                  j4::AngularType, j5::AngularType, j6::AngularType)
    dj1, dj2, dj3, dj4, dj5, dj6 = Int64.((2j1, 2j2, 2j3, 2j4, 2j5, 2j6))
    check_couple(dj1, dj2, dj3) & check_couple(dj1, dj5, dj6) &
    check_couple(dj4, dj2, dj6) & check_couple(dj4, dj5, dj3)
end

# 测试 9j 系数的特殊情况
@testset "9j" begin
for j1 = test_collect2, j2 = test_collect2, j3 = test_collect2,
    j4 = test_collect2, j5 = test_collect2, j7 = test_collect2
    if check_6j(j1, j2, j3, j5, j4, j7)
        nj = nineJ(j1, j2, j3, j4, j5, j3, j7, j7, 0)
        snj = iphase(Int(j2+j3+j4+j7)) / SqrtRational((2j3+1)*(2j7+1))
        snj *= sixJ(j1, j2, j3, j5, j4, j7)
        @test nj == snj
    end
end
end