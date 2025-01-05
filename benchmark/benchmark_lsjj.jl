using CGcoefficient
using BenchmarkTools

function calculate_lsjj(test_range::AbstractArray)
    sum = 0.0
    for l1 in test_range
        for l2 in test_range
            for (dj1, dj2) in [(2l1-1, 2l2-1), (2l1-1, 2l2+1), (2l1+1, 2l2-1), (2l1+1, 2l2+1)]
                for L in abs(l1-l2):(l1+l2)
                    for (S, J) in [(0, L), (1, L-1), (1, L), (1, L+1)]
                        sum += flsjj(l1, l2, dj1, dj2, L, S, J)
                    end
                end
            end
        end
    end
    return sum
end

function calculate_norm9j(test_range::AbstractArray)
    sum = 0.0
    for l1 in test_range
        for l2 in test_range
            for (dj1, dj2) in [(2l1-1, 2l2-1), (2l1-1, 2l2+1), (2l1+1, 2l2-1), (2l1+1, 2l2+1)]
                for L in abs(l1-l2):(l1+l2)
                    for (S, J) in [(0, L), (1, L-1), (1, L), (1, L+1)]
                        sum += fnorm9j(2l1, 1, dj1, 2l2, 1, dj2, 2L, 2S, 2J)
                    end
                end
            end
        end
    end
    return sum
end

test_range = 0:10

wigner_init_float(maximum(test_range), "Jmax", 9)

t1 = @belapsed calculate_lsjj(test_range)
t2 = @belapsed calculate_norm9j(test_range)

println("diff = ", calculate_lsjj(test_range) - calculate_norm9j(test_range))
println("lsjj time = $(t1*1000)ms")
println("norm9j time = $(t2*1000)ms")