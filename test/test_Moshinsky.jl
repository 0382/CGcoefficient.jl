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

const inv_sqrt2 = exact_sqrt(1//2)
const Moshinsky_test_set_result = [
    1, inv_sqrt2, -inv_sqrt2, 1//2, -inv_sqrt2, 1//2, inv_sqrt2, 0, -inv_sqrt2, -1, inv_sqrt2, 0, -inv_sqrt2 
]

function test_Moshinsky()
    wigner_init_float(10, "nmax", 0)
    for i = eachindex(Moshinsky_test_set_result)
        N, L, n, l, n1, l1, n2, l2, Λ = Moshinsky_test_set[i, :]
        @test Moshinsky(N, L, n, l, n1, l1, n2, l2, Λ) == Moshinsky_test_set_result[i]
        @test fMoshinsky(N, L, n, l, n1, l1, n2, l2, Λ) ≈ float(Moshinsky_test_set_result[i])
    end
end