const _fbinomial_data = Ref(Float64[])
const _fbinomial_nmax = Ref{Int}(0)

@inline get_fbinomial_data() = _fbinomial_data[]
@inline get_fbinomial_nmax() = _fbinomial_nmax[]

@inline function fbinomial_data_size(n::Int)
    x = div(n, 2) + 1
    return x * (x + isodd(n))
end

@inline function fbinomial_index(n::Int, k::Int)
    x = div(n, 2) + 1
    return x * (x - iseven(n)) + k + 1
end

"""
    fbinomial(n::Integer, k::Integer)
`binomial` with Float64 return value.
"""
function fbinomial(n::Integer, k::Integer)::Float64
    if n < 0 || n > get_fbinomial_nmax() || k < 0 || k > n
        return 0.0
    end
    if k > div(n, 2)
        k = n - k
    end
    return @inbounds get_fbinomial_data()[fbinomial_index(n, k)]
end

function _extent_fbinomial_data(n::Int)
    if n > get_fbinomial_nmax()
        old_data = copy(get_fbinomial_data())
        global _fbinomial_data[] = Vector{Float64}(undef, fbinomial_data_size(n))
        copyto!(get_fbinomial_data(), old_data)
        for m = get_fbinomial_nmax()+1:n
            for k = 0:div(m, 2)
                get_fbinomial_data()[fbinomial_index(m, k)] = fbinomial(m - 1, k) + fbinomial(m - 1, k - 1)
            end
            global _fbinomial_nmax[] = get_fbinomial_nmax() + 1
        end
        global _fbinomial_nmax[] = n
    end
    nothing
end

"""
    reserve_fbinomial(n::Integer, mode::AbstractString, rank::Integer)
This function reserves memory for fbinomial(n, k). In this code, the `fbinomial` function is only valid in the stored range. If you call a `fbinomial` function out of the range, it just gives you `0`. The `__init__()` function stores to `nmax = 67`.

The parameters means

|                                       |    Calculate range    |   CG & 3j   | 6j & Racah  |     9j      |
| :-----------------------------------: | :-------------------: | :---------: | :---------: | :---------: |
|          meaning of `type`            | `type`\\\\`rank`      |      3      |      6      |      9      |
|        max angular momentum           |    `"Jmax"`           | `3*Jmax+1`  | `4*Jmax+1`  | `5*Jmax+1`  |
| max two-body coupled angular momentum |   `"2bjmax"`          | `2*jmax+1`  | `3*jmax+1`  | `4*jmax+1`  |
|            max binomial               |    `"nmax"`           |   `nmax`    |   `namx`    |   `nmax`    |

The `"2bjmax"` mode means your calculation only consider two-body coupling, and no three-body coupling. This mode assumes that in all these coefficients, at least one of the angular momentun is just a single particle angular momentum. With this assumption, `"2bjmax"` mode will use less memory than `"Jmax"` mode.

`"Jmax"` means the global maximum angular momentum, for every parameters. It is always safe with out any assumption.

The `"nmax"` mode directly set `nmax`, and the `rank` parameter is ignored. 

`rank = 6` means you only need to calculate CG and/or 6j symbols, you don't need to calculate 9j symbol.

For example
```julia
reserve_fbinomial(21, "Jmax", 6)
```
means you calculate CG and 6j symbols, and donot calculate 9j symbol. The maximum angular momentum in your system is 21.

You do not need to rememmber those values in the table. You just need to find the maximum angular momentum in you canculation, then call the function.

The reserve_fbinomial function is **not** thread safe, so you should call it before you start your calculation.
"""
function reserve_fbinomial(n::Integer, mode::AbstractString, rank::Integer)
    if mode == "Jmax"
        rank == 3 && _extent_fbinomial_data(3 * n + 1)
        rank == 6 && _extent_fbinomial_data(4 * n + 1)
        rank == 9 && _extent_fbinomial_data(5 * n + 1)
        if rank != 3 && rank != 6 && rank != 9
            throw(ArgumentError("invalid rank $rank"))
        end
    elseif mode == "2bjmax"
        rank == 3 && _extent_fbinomial_data(2 * n + 1)
        rank == 6 && _extent_fbinomial_data(3 * n + 1)
        rank == 9 && _extent_fbinomial_data(4 * n + 1)
        if rank != 3 && rank != 6 && rank != 9
            throw(ArgumentError("invalid rank $rank"))
        end
    elseif mode == "nmax"
        _extent_fbinomial_data(n)
    else
        throw(ArgumentError("invalid mode $mode"))
    end
    nothing
end


# basic CG coefficient calculation function with Float64 return value
# unlike the SqrtRational version, this function gives out zero when the parameters are not valid
function _fCG(dj1::Int64, dj2::Int64, dj3::Int64, dm1::Int64, dm2::Int64, dm3::Int64)
    check_jm(dj1, dm1) && check_jm(dj2, dm2) && check_jm(dj3, dm3) || return zero(Float64)
    check_couple(dj1, dj2, dj3) || return zero(Float64)
    dm1 + dm2 == dm3 || return zero(Float64)
    J::Int64 = div(dj1 + dj2 + dj3, 2)
    Jm1::Int64 = J - dj1
    Jm2::Int64 = J - dj2
    Jm3::Int64 = J - dj3
    j1mm1::Int64 = div(dj1 - dm1, 2)
    j2mm2::Int64 = div(dj2 - dm2, 2)
    j3mm3::Int64 = div(dj3 - dm3, 2)
    j2pm2::Int64 = div(dj2 + dm2, 2)
    A = sqrt(fbinomial(dj1, Jm2) * fbinomial(dj2, Jm3) / (
        fbinomial(J + 1, Jm3) * fbinomial(dj1, j1mm1) *
        fbinomial(dj2, j2mm2) * fbinomial(dj3, j3mm3)
    ))
    B = zero(Float64)
    low::Int64 = max(zero(Int64), j1mm1 - Jm2, j2pm2 - Jm1)
    high::Int64 = min(Jm3, j1mm1, j2pm2)
    for z in low:high
        B = -B + fbinomial(Jm3, z) * fbinomial(Jm2, j1mm1 - z) * binomial(Jm1, j2pm2 - z)
    end
    return iphase(high) * A * B
end

function _f6j(dj1::Int64, dj2::Int64, dj3::Int64, dj4::Int64, dj5::Int64, dj6::Int64)
    check_couple(dj1, dj2, dj3) && check_couple(dj1, dj5, dj6) & check_couple(dj4, dj2, dj6) && check_couple(dj4, dj5, dj3) || return zero(Float64)
    j123::Int64 = div(dj1 + dj2 + dj3, 2)
    j156::Int64 = div(dj1 + dj5 + dj6, 2)
    j426::Int64 = div(dj4 + dj2 + dj6, 2)
    j453::Int64 = div(dj4 + dj5 + dj3, 2)
    jpm123::Int64 = div(dj1 + dj2 - dj3, 2)
    jpm132::Int64 = div(dj1 + dj3 - dj2, 2)
    jpm231::Int64 = div(dj2 + dj3 - dj1, 2)
    jpm156::Int64 = div(dj1 + dj5 - dj6, 2)
    jpm426::Int64 = div(dj4 + dj2 - dj6, 2)
    jpm453::Int64 = div(dj4 + dj5 - dj3, 2)
    A = sqrt(fbinomial(j123 + 1, dj1 + 1) * fbinomial(dj1, jpm123) / (
        fbinomial(j156 + 1, dj1 + 1) * fbinomial(dj1, jpm156) *
        fbinomial(j453 + 1, dj4 + 1) * fbinomial(dj4, jpm453) *
        fbinomial(j426 + 1, dj4 + 1) * fbinomial(dj4, jpm426)
    ))
    B = zero(Float64)
    low::Int64 = max(j123, j453, j426, j156)
    high::Int64 = min(jpm123 + j453, jpm132 + j426, jpm231 + j156)
    for x = low:high
        B = -B + fbinomial(x + 1, j123 + 1) * fbinomial(jpm123, x - j453) *
                 fbinomial(jpm132, x - j426) * fbinomial(jpm231, x - j156)
    end
    return iphase(high) * A * B / (dj4 + 1)
end

function _f9j(dj1::Int64, dj2::Int64, dj3::Int64,
    dj4::Int64, dj5::Int64, dj6::Int64,
    dj7::Int64, dj8::Int64, dj9::Int64)
    check_couple(dj1, dj2, dj3) && check_couple(dj4, dj5, dj6) && check_couple(dj7, dj8, dj9) || return zero(Float64)
    check_couple(dj1, dj4, dj7) && check_couple(dj2, dj5, dj8) && check_couple(dj3, dj6, dj9) || return zero(Float64)
    j123::Int64 = div(dj1 + dj2 + dj3, 2)
    j456::Int64 = div(dj4 + dj5 + dj6, 2)
    j789::Int64 = div(dj7 + dj8 + dj9, 2)
    j147::Int64 = div(dj1 + dj4 + dj7, 2)
    j258::Int64 = div(dj2 + dj5 + dj8, 2)
    j369::Int64 = div(dj3 + dj6 + dj9, 2)
    pm123::Int64 = div(dj1 + dj2 - dj3, 2)
    pm132::Int64 = div(dj1 + dj3 - dj2, 2)
    pm231::Int64 = div(dj2 + dj3 - dj1, 2)
    pm456::Int64 = div(dj4 + dj5 - dj6, 2)
    pm465::Int64 = div(dj4 + dj6 - dj5, 2)
    pm564::Int64 = div(dj5 + dj6 - dj4, 2)
    pm789::Int64 = div(dj7 + dj8 - dj9, 2)
    pm798::Int64 = div(dj7 + dj9 - dj8, 2)
    pm897::Int64 = div(dj8 + dj9 - dj7, 2)
    # value in sqrt
    P0_nu = fbinomial(j123 + 1, dj1 + 1) * fbinomial(dj1, pm123) *
            fbinomial(j456 + 1, dj5 + 1) * fbinomial(dj5, pm456) *
            fbinomial(j789 + 1, dj9 + 1) * fbinomial(dj9, pm798)
    P0_de = fbinomial(j147 + 1, dj1 + 1) * fbinomial(dj1, div(dj1 + dj4 - dj7, 2)) *
            fbinomial(j258 + 1, dj5 + 1) * fbinomial(dj5, div(dj2 + dj5 - dj8, 2)) *
            fbinomial(j369 + 1, dj9 + 1) * fbinomial(dj9, div(dj3 + dj9 - dj6, 2))
    P0 = P0_nu / P0_de
    dtl::Int64 = max(abs(dj2 - dj6), abs(dj4 - dj8), abs(dj1 - dj9))
    dth::Int64 = min(dj2 + dj6, dj4 + dj8, dj1 + dj9)
    PABC = zero(Float64)
    for dt::Int64 = dtl:2:dth
        j19t::Int64 = div(dj1 + dj9 + dt, 2)
        j26t::Int64 = div(dj2 + dj6 + dt, 2)
        j48t::Int64 = div(dj4 + dj8 + dt, 2)
        Pt_de = fbinomial(j19t + 1, dt + 1) * fbinomial(dt, div(dj1 + dt - dj9, 2)) *
                fbinomial(j26t + 1, dt + 1) * fbinomial(dt, div(dj2 + dt - dj6, 2)) *
                fbinomial(j48t + 1, dt + 1) * fbinomial(dt, div(dj4 + dt - dj8, 2))
        Pt_de *= (dt + 1)^2
        xl::Int64 = max(j123, j369, j26t, j19t)
        xh::Int64 = min(pm123 + j369, pm132 + j26t, pm231 + j19t)
        At = zero(Float64)
        for x = xl:xh
            At = -At + fbinomial(x + 1, j123 + 1) * fbinomial(pm123, x - j369) *
                       fbinomial(pm132, x - j26t) * fbinomial(pm231, x - j19t)
        end
        yl::Int64 = max(j456, j26t, j258, j48t)
        yh::Int64 = min(pm456 + j26t, pm465 + j258, pm564 + j48t)
        Bt = zero(Float64)
        for y = yl:yh
            Bt = -Bt + fbinomial(y + 1, j456 + 1) * fbinomial(pm456, y - j26t) *
                       fbinomial(pm465, y - j258) * fbinomial(pm564, y - j48t)
        end
        zl::Int64 = max(j789, j19t, j48t, j147)
        zh::Int64 = min(pm789 + j19t, pm798 + j48t, pm897 + j147)
        Ct = zero(Float64)
        for z = zl:zh
            Ct = -Ct + fbinomial(z + 1, j789 + 1) * fbinomial(pm789, z - j19t) *
                       fbinomial(pm798, z - j48t) * fbinomial(pm897, z - j147)
        end
        PABC += (iphase(xh + yh + zh) * At * Bt * Ct) / Pt_de
    end
    return iphase(dth) * sqrt(P0) * PABC
end

"""
    fCG(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer)
float64 and fast CG coefficient.
"""
function fCG(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer)
    return _fCG(Int64.((dj1, dj2, dj3, dm1, dm2, dm3))...)
end

"""
    f3j(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer)
float64 and fast Wigner 3j symbol.
"""
function f3j(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer)
    return iphase(dj1 + div(dj3 + dm3, 2)) * fCG(dj1, dj2, dj3, -dm1, -dm2, dm3) / sqrt(dj3 + 1)
end

"""
    f6j(dj1::Integer, dj2::Integer, dj3::Integer, dj4::Integer, dj5::Integer, dj6::Integer)
float64 and fast Wigner 6j symbol.
"""
function f6j(dj1::Integer, dj2::Integer, dj3::Integer, dj4::Integer, dj5::Integer, dj6::Integer)
    return _f6j(Int64.((dj1, dj2, dj3, dj4, dj5, dj6))...)
end

"""
    Racah(dj1::Integer, dj2::Integer, dj3::Integer, dj4::Integer, dj5::Integer, dj6::Integer)
float64 and fast Racah coefficient.
"""
function fRacah(dj1::Integer, dj2::Integer, dj3::Integer, dj4::Integer, dj5::Integer, dj6::Integer)
    return iphase(div(dj1 + dj2 + dj3 + dj4, 2)) * f6j(dj1, dj2, dj5, dj4, dj3, dj6)
end

"""
    f9j(dj1::Integer, dj2::Integer, dj3::Integer,
        dj4::Integer, dj5::Integer, dj6::Integer,
        dj7::Integer, dj8::Integer, dj9::Integer)
float64 and fast Wigner 9j symbol.
"""
function f9j(dj1::Integer, dj2::Integer, dj3::Integer,
             dj4::Integer, dj5::Integer, dj6::Integer,
             dj7::Integer, dj8::Integer, dj9::Integer)
    return _f9j(Int64.((dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9))...)
end
