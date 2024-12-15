# benchmark the exact version and the float version

using CGcoefficient
using BenchmarkTools
using WignerSymbols

function bench_CG(func::Function, djmax::Int)
    sum = 0.0
    for dj1 in 1:djmax
        print(dj1, '\r')
        for dj2 in 0:dj1
            for dj3 = abs(dj1-dj2):2:(dj1+dj2)
                for dm1 in -dj1:2:0
                    for dm2 in -dj2:2:dj2
                        dm3 = dm1+dm2
                        ans = func(dj1, dj2, dj3, dm1, dm2, dm3)
                        sum += ans
                    end
                end
            end 
        end
    end
    return sum
end

wigner_init_exact(100, "Jmax", 3)
wigner_init_float(100, "Jmax", 3)
