using CGcoefficient
using WignerSymbols

function get_max_error_3j(djmax)
    wigner_init_float(div(djmax, 2) + 1, "Jmax", 3)
    max_error = 0.0
    max_rel_error = 0.0
    for dj1 in (djmax & 1):2:djmax
    for dj2 in dj1:2:djmax
    for dj3 in (dj2-dj1):2:djmax
        if dj1 < djmax && dj2 < djmax && dj3 < (djmax-isodd(djmax)) # bacause they are calculated in privious call
            continue
        end
        for dm1 in -dj1:2:dj1
            for dm2 in (dj2 & 1):2:dj2 # minus dm2 has same absolute value of (-dm1, -dm2)
                dm3 = -dm1-dm2
                check_jm(dj3, dm3) || continue
                exact = wigner3j(dj1//2, dj2//2, dj3//2, dm1//2, dm2//2)
                fval = f3j(dj1, dj2, dj3, dm1, dm2, dm3)
                max_error = max(max_error, abs(exact - fval))
                if exact != 0.0
                    max_rel_error = max(max_rel_error, abs((exact - fval)/exact))
                end
            end
        end
    end end end
    return max_error, max_rel_error
end

err = zeros(100)
rel_err = zeros(100)

for dj in 1:100
    err[dj], rel_err[dj] = @time get_max_error_3j(dj)
    println("dj = $dj, err = $(err[dj]), rel_err = $(rel_err[dj])")
end