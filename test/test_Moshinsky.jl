# this file test the Moshinsky coefficients
# the data comes from: Buck et al. Nuc. Phys. A 600 (1996) 387-402

const Moshinsky_test_set = [
    0 0 0 0 0 0 0 0 0
    0 1 0 0 0 0 0 1 1
    0 0 0 1 0 0 0 1 1
    0 2 0 0 0 0 0 2 2
    0 1 0 1 0 0 0 2 2
    0 0 0 2 0 0 0 2 2
    1 0 0 0 0 1 0 1 0
    0 1 0 1 0 1 0 1 0
    0 0 1 0 0 1 0 1 0
    0 1 0 1 0 1 0 1 1
    0 2 0 0 0 1 0 1 2
    0 1 0 1 0 1 0 1 2
    0 0 0 2 0 1 0 1 2]

const Moshinsky_test_set2 = [
    0 2 1 0 0 1 0 3 2
    0 1 0 5 0 1 0 5 6
    0 1 0 3 0 2 0 2 4
    1 3 0 1 0 2 0 4 3
    0 5 0 2 0 2 0 5 4
    0 3 1 6 2 2 1 3 4
    1 0 2 5 2 2 1 3 5
    0 2 4 2 2 2 1 4 2
    3 2 0 4 2 2 1 4 4
]

const Moshinsky_test_set_result2 = [
    SqrtRational(-1//2, 7//10),
    SqrtRational(1//2, 1),
    SqrtRational(0, 1),
    SqrtRational(-1//2, 5//14),
    SqrtRational(-1//6, 1//2),
    SqrtRational(1//24, 65//3),
    SqrtRational(-5//96, 13//7),
    SqrtRational(-1//28, 195//14),
    SqrtRational(4463//25872, 1//3)
]

const inv_sqrt2 = exact_sqrt(1//2)
const Moshinsky_test_set_result = [
    1, inv_sqrt2, -inv_sqrt2, 1//2, -inv_sqrt2, 1//2, inv_sqrt2, 0, -inv_sqrt2, -1, inv_sqrt2, 0, -inv_sqrt2 
]

function test_Moshinsky()
    wigner_init_float(7, "Moshinsky", 0)
    for i = eachindex(Moshinsky_test_set_result)
        N, L, n, l, n1, l1, n2, l2, Λ = Moshinsky_test_set[i, :]
        @test Moshinsky(N, L, n, l, n1, l1, n2, l2, Λ) == Moshinsky_test_set_result[i]
        @test fMoshinsky(N, L, n, l, n1, l1, n2, l2, Λ) ≈ float(Moshinsky_test_set_result[i])
    end
    for i = eachindex(Moshinsky_test_set_result2)
        N, L, n, l, n1, l1, n2, l2, Λ = Moshinsky_test_set2[i, :]
        @test Moshinsky(N, L, n, l, n1, l1, n2, l2, Λ) == Moshinsky_test_set_result2[i]
        @test fMoshinsky(N, L, n, l, n1, l1, n2, l2, Λ) ≈ float(Moshinsky_test_set_result2[i])
    end
end