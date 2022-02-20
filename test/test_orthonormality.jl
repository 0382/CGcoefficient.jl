# this file test the orthonormality of Wigner symbols
# All formulat are refered in
# [1] D. A. Varshalovich, A. N. Moskalev and V. K. Khersonskii, Quantum Theory of Angular Momentum, (World Scientific, 1988).

# Ref[1], P236, Sec 8.1, Formula (8)
function test_orthonormality_CG(test_range::AbstractArray)
    # Formula (8) first part
    for j1 in test_range
    for j2 in test_range
    for j3 in test_range
    for j4 in test_range
        if check_couple(Int.((2j1, 2j2, 2j3))...) && check_couple(Int.((2j1, 2j2, 2j4))...)
            for m3 = -min(j3,j4):1:min(j3,j4)
                sum = zero(SqrtRational)
                for m1 = -j1:1:j1
                    m2 = m3 - m1
                    if check_jm(Int(2j2), Int(2m2))
                        sum += CG(j1,j2,j3,m1,m2,m3) * CG(j1,j2,j4,m1,m2,m3)
                    end
                end
                @test sum == Int(j3==j4)
            end
        end
    end end end end
    # Formula (8) second part
    for j1 in test_range
    for j2 in test_range
        for m1 = -j1:1:j1
        for m2 = -j2:1:j2
        for m3 = -j1:1:j1
            m = m1 + m2
            m4 = m - m3
            sum = zero(SqrtRational)
            for j3 = abs(j1-j2):1:(j1+j2)
                if check_jm(Int(2j2), Int(2m4)) && check_jm(Int(2j3), Int(2m))
                    sum += CG(j1,j2,j3,m1,m2,m) * CG(j1,j2,j3,m3,m4,m)
                end
            end
            @test sum == Int(m1 == m3)
        end end end
    end end
end

# Ref[1], P453, Sec 12.1, Formula (2)
function test_summation_3j(test_range::AbstractArray)
    for j1 = test_range
    for j3 = 0:1//2:maximum(test_range)
        check_couple(Int(2j1), Int(2j1), Int(2j3)) || continue
        sum = zero(SqrtRational)
        for m1 = -j1:1:j1
            sum += iphase(Int(j1-m1)) * threeJ(j1,j1,j3,m1,-m1,0)
        end
        @test sum == exact_sqrt(2j1+1) * Int(j3==0)
    end end
end

# Ref[1], P291, Sec 9.1, Formula(9)
function test_orthonormality_6j(test_range::AbstractArray)
    for j1 in test_range
    for j2 in test_range
    for j4 in test_range
    for j5 in test_range
    for j6 in test_range
    for j7 in test_range
        is_same_parity(Int(2j6), Int(2j7)) || continue
        check_couple(Int(2j1), Int(2j5), Int(2j6)) || continue
        check_couple(Int(2j1), Int(2j5), Int(2j7)) || continue
        check_couple(Int(2j4), Int(2j2), Int(2j6)) || continue
        check_couple(Int(2j4), Int(2j2), Int(2j7)) || continue
        is_same_parity(Int(2(j1+j2)), Int(2(j4+j5))) || continue
        sum = zero(SqrtRational)
        jmin = max(abs(j1-j2), abs(j4-j5))
        jmax = min(j1+j2, j4+j5)
        for j3 = jmin:1:jmax
            sum += (2j3+1) * (2j6+1) * (
                sixJ(j1,j2,j3,j4,j5,j6) * sixJ(j1,j2,j3,j4,j5,j7)
            )
        end
        @test sum == Int(j6==j7)
    end end end end end end
end

# Ref[1], P335, Sec 10.1, Formula (9)
function test_orthonormality_9j(test_range::AbstractArray)
    for a in test_range
    for b in test_range
    for c in test_range
    for d in test_range
    for e in test_range
    for f in test_range
    for cx in test_range
    for fx in test_range
    for j in test_range
        is_same_parity(Int(2c), Int(2cx)) || continue
        is_same_parity(Int(2f), Int(2fx)) || continue
        check_couple(Int(2a), Int(2b), Int(2c)) || continue
        check_couple(Int(2a), Int(2b), Int(2cx)) || continue
        check_couple(Int(2d), Int(2e), Int(2f)) || continue
        check_couple(Int(2d), Int(2e), Int(2fx)) || continue
        check_couple(Int(2c), Int(2f), Int(2j)) || continue
        check_couple(Int(2cx), Int(2fx), Int(2j)) || continue
        sum = zero(SqrtRational)
        for g in abs(a-d):1:(a+d)
        for h in abs(b-e):1:(b+e)
            if check_couple(Int(2g), Int(2h), Int(2j))
                sum += (2g+1)*(2h+1)*(
                    nineJ(a,b,c,d,e,f,g,h,j) * nineJ(a,b,cx,d,e,fx,g,h,j)
                )
            end
        end end
        if (c==cx) && (f==fx)
            @test sum == 1//((2*c+1)*(2*f+1))
        else
            @test sum == 0
        end
    end end end end end end end end end
end