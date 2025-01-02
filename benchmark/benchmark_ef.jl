# benchmark the exact version and the float version

using CGcoefficient

function bench_3j(func::Function, djmax::Int)
    sum = 0.0
    for dj1 in 1:djmax
        print(dj1, '\r')
        for dj2 in 0:dj1
            for dj3 = abs(dj1-dj2):2:dj2
                for dm1 in -dj1:2:0
                    for dm2 in -dj2:2:dj2
                        dm3 = -dm1-dm2
                        ans = func(dj1, dj2, dj3, dm1, dm2, dm3)
                        sum += ans
                    end
                end
            end 
        end
    end
    return sum
end

function bench_6j(func::Function, djmax::Int)
    sum = 0.0
    for dj1 in 0:djmax
        print(dj1, '\r')
        for dj2 in 0:dj1
            for dj3 in abs(dj1-dj2):2:dj2
                for dj4 in 0:dj1
                    for dj5 in abs(dj3-dj4):2:min(dj1, dj3+dj4)
                        dj6min = max(abs(dj1 - dj5), abs(dj2 - dj4))
                        dj6max = min(dj1, dj2 + dj4)
                        for dj6 in dj6min:2:dj6max
                            ans = func(dj1, dj2, dj3, dj4, dj5, dj6)
                            sum += ans
                        end
                    end
                end
            end
        end
    end
    return sum
end

Jmax = 20

wigner_init_pf(Jmax, "Jmax", 6)
wigner_init_float(Jmax, "Jmax", 6)

big3j(dj1,dj2,dj3,dm1,dm2,dm3) = convert(Float64, float(d3j(dj1,dj2,dj3,dm1,dm2,dm3)))
big6j(dj1,dj2,dj3,dj4,dj5,dj6) = convert(Float64, float(d6j(dj1,dj2,dj3,dj4,dj5,dj6)))

bench_3j(f3j, 1)
bench_3j(ef3j, 1)
bench_3j(big3j, 1)
bench_6j(f6j, 1)
bench_6j(ef6j, 1)
bench_6j(big6j, 1)

t_f3j = @elapsed bench_3j(f3j, 2Jmax)
t_ef3j = @elapsed bench_3j(ef3j, 2Jmax)
t_big3j = @elapsed bench_3j(big3j, 2Jmax)

t_f6j = @elapsed bench_6j(f6j, 2Jmax)
t_ef6j = @elapsed bench_6j(ef6j, 2Jmax)
t_big6j = @elapsed bench_6j(big6j, 2Jmax)

println("f3j: ", t_f3j)
println("ef3j: ", t_ef3j)
println("big3j: ", t_big3j)
println("f6j: ", t_f6j)
println("ef6j: ", t_ef6j)
println("big6j: ", t_big6j)
