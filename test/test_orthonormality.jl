# this file test the orthonormality of Wigner symbols
# All formulat are refered in
# [1] D. A. Varshalovich, A. N. Moskalev and V. K. Khersonskii, Quantum Theory of Angular Momentum, (World Scientific, 1988).

# Ref[1], P236, Sec 8.1, Formula (8)
function test_orthonormality_dCG(test_range::AbstractArray)
    # Formula (8) first part
    for dj1 in test_range
    for dj2 in test_range
    for dj3 in test_range
    for dj4 in test_range
        if check_couple(dj1, dj2, dj3) && check_couple(dj1, dj2, dj4)
            for dm3 = -min(dj3,dj4):2:min(dj3,dj4)
                sum = zero(SqrtRational)
                for dm1 = -dj1:2:dj1
                    dm2 = dm3 - dm1
                    if check_jm(dj2, dm2)
                        sum += dCG(dj1,dj2,dj3,dm1,dm2,dm3) * dCG(dj1,dj2,dj4,dm1,dm2,dm3)
                    end
                end
                @test sum == Int(dj3==dj4)
            end
        end
    end end end end
    # Formula (8) second part
    for dj1 in test_range
    for dj2 in test_range
        for dm1 = -dj1:2:dj1
        for dm2 = -dj2:2:dj2
        for dm3 = -dj1:2:dj1
            dm = dm1 + dm2
            dm4 = dm - dm3
            sum = zero(SqrtRational)
            for dj3 = abs(dj1-dj2):2:(dj1+dj2)
                if check_jm(dj2, dm4) && check_jm(dj3, dm)
                    sum += dCG(dj1,dj2,dj3,dm1,dm2,dm) * dCG(dj1,dj2,dj3,dm3,dm4,dm)
                end
            end
            @test sum == Int(dm1 == dm3)
        end end end
    end end
end

# Ref[1], P453, Sec 12.1, Formula (2)
function test_summation_d3j(test_range::AbstractArray)
    for dj1 = test_range
    for dj3 = 0:maximum(test_range)
        check_couple(dj1, dj1, dj3) || continue
        sum = zero(SqrtRational)
        for dm1 = -dj1:2:dj1
            sum += iphase(div(dj1-dm1, 2)) * d3j(dj1,dj1,dj3,dm1,-dm1,0)
        end
        @test sum == exact_sqrt(dj1+1) * Int(dj3==0)
    end end
end

# Ref[1], P291, Sec 9.1, Formula(9)
function test_orthonormality_d6j(test_range::AbstractArray)
    for dj1 in test_range
    for dj2 in test_range
    for dj4 in test_range
    for dj5 in test_range
    for dj6 in test_range
    for dj7 in test_range
        is_same_parity(dj6, dj7) || continue
        check_couple(dj1, dj5, dj6) || continue
        check_couple(dj1, dj5, dj7) || continue
        check_couple(dj4, dj2, dj6) || continue
        check_couple(dj4, dj2, dj7) || continue
        is_same_parity(dj1 + dj2, dj4 + dj5) || continue
        sum = zero(SqrtRational)
        djmin = max(abs(dj1-dj2), abs(dj4-dj5))
        djmax = min(dj1+dj2, dj4+dj5)
        for dj3 = djmin:2:djmax
            sum += (dj3+1) * (dj6+1) * (
                d6j(dj1,dj2,dj3,dj4,dj5,dj6) * d6j(dj1,dj2,dj3,dj4,dj5,dj7)
            )
        end
        @test sum == Int(dj6==dj7)
    end end end end end end
end

# Ref[1], P335, Sec 10.1, Formula (9)
function test_orthonormality_d9j(test_range::AbstractArray)
    for a in test_range
    for b in test_range
    for c in test_range
    for d in test_range
    for e in test_range
    for f in test_range
    for cx in test_range
    for fx in test_range
    for j in test_range
        is_same_parity(c, cx) || continue
        is_same_parity(f, fx) || continue
        check_couple(a, b, c) || continue
        check_couple(a, b, cx) || continue
        check_couple(d, e, f) || continue
        check_couple(d, e, fx) || continue
        check_couple(c, f, j) || continue
        check_couple(cx, fx, j) || continue
        sum = zero(SqrtRational)
        for g in abs(a-d):2:(a+d)
        for h in abs(b-e):2:(b+e)
            if check_couple(g, h, j)
                sum += (g+1)*(h+1)*(
                    d9j(a,b,c,d,e,f,g,h,j) * d9j(a,b,cx,d,e,fx,g,h,j)
                )
            end
        end end
        if (c==cx) && (f==fx)
            @test sum == 1//((c+1)*(f+1))
        else
            @test sum == 0
        end
    end end end end end end end end end
end