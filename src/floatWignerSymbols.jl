_fbinomial_data = Float64[]
_fbinomial_nmax = Int(0)
@inline get_fbinomial_data()::Vector{Float64} = _fbinomial_data::Vector{Float64}
@inline get_fbinomial_nmax()::Int = _fbinomial_nmax::Int

"""
    fbinomial(n::Integer, k::Integer)
`binomial` with Float64 return value.
"""
@inline function fbinomial(n::Integer, k::Integer)::Float64
    if n < 0 || n > get_fbinomial_nmax() || k < 0 || k > n
        return 0.0
    end
    k = min(k, n - k)
    return @inbounds get_fbinomial_data()[binomial_index(n, k)]
end

"""
    function unsafe_fbinomial(n::Int, k::Int)
Same as fbinomial, but without bounds check. Thus for `n < 0` or `n > _nmax` or `k < 0` or `k > n`, the result is undefined.
In Wigner symbol calculation, the mathematics guarantee that `unsafe_fbinomial(n, k)` is always safe.
"""
@inline function unsafe_fbinomial(n::Int, k::Int)::Float64
    k = min(k, n - k)
    return @inbounds get_fbinomial_data()[binomial_index(n, k)]
end

function _extent_fbinomial_data(n::Int)
    if n > get_fbinomial_nmax()
        old_data = copy(get_fbinomial_data())
        resize!(get_fbinomial_data(), binomial_data_size(n))
        copyto!(get_fbinomial_data(), old_data)
        for m = get_fbinomial_nmax()+1:n
            for k = 0:div(m, 2)
                get_fbinomial_data()[binomial_index(m, k)] = binomial(BigInt(m), BigInt(k))
            end
            global _fbinomial_nmax = get_fbinomial_nmax() + 1
        end
        global _fbinomial_nmax = n
    end
    nothing
end

"""
    wigner_init_float(n::Integer, mode::AbstractString, rank::Integer)
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
wigner_init_float(21, "Jmax", 6)
```
means you calculate CG and 6j symbols, and donot calculate 9j symbol. The maximum angular momentum in your system is 21.

You do not need to rememmber those values in the table. You just need to find the maximum angular momentum in you canculation, then call the function.

The wigner_init_float function is **not** thread safe, so you should call it before you start your calculation.
"""
function wigner_init_float(n::Integer, mode::AbstractString, rank::Integer)
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
    A = sqrt(unsafe_fbinomial(dj1, Jm2) * unsafe_fbinomial(dj2, Jm3) / (
        unsafe_fbinomial(J + 1, Jm3) * unsafe_fbinomial(dj1, j1mm1) *
        unsafe_fbinomial(dj2, j2mm2) * unsafe_fbinomial(dj3, j3mm3)
    ))
    B = zero(Float64)
    low::Int64 = max(zero(Int64), j1mm1 - Jm2, j2pm2 - Jm1)
    high::Int64 = min(Jm3, j1mm1, j2pm2)
    for z in low:high
        B = -B + unsafe_fbinomial(Jm3, z) * unsafe_fbinomial(Jm2, j1mm1 - z) * binomial(Jm1, j2pm2 - z)
    end
    return iphase(high) * A * B
end

function _fCG0(j1::Int64, j2::Int64, j3::Int64)
    check_couple(2j1, 2j2, 2j3) || return zero(Float64)
    J = j1 + j2 + j3
    isodd(J) && return zero(Float64)
    g = div(J, 2)
    return iphase(g - j3) * unsafe_fbinomial(g, j3) * unsafe_fbinomial(j3, g - j1) / sqrt(unsafe_fbinomial(J + 1, 2j3 + 1) * unsafe_fbinomial(2j3, J - 2j1))
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
    A = sqrt(unsafe_fbinomial(j123 + 1, dj1 + 1) * unsafe_fbinomial(dj1, jpm123) / (
        unsafe_fbinomial(j156 + 1, dj1 + 1) * unsafe_fbinomial(dj1, jpm156) *
        unsafe_fbinomial(j453 + 1, dj4 + 1) * unsafe_fbinomial(dj4, jpm453) *
        unsafe_fbinomial(j426 + 1, dj4 + 1) * unsafe_fbinomial(dj4, jpm426)
    ))
    B = zero(Float64)
    low::Int64 = max(j123, j453, j426, j156)
    high::Int64 = min(jpm123 + j453, jpm132 + j426, jpm231 + j156)
    for x = low:high
        B = -B + unsafe_fbinomial(x + 1, j123 + 1) * unsafe_fbinomial(jpm123, x - j453) *
                 unsafe_fbinomial(jpm132, x - j426) * unsafe_fbinomial(jpm231, x - j156)
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
    P0_nu = unsafe_fbinomial(j123 + 1, dj1 + 1) * unsafe_fbinomial(dj1, pm123) *
            unsafe_fbinomial(j456 + 1, dj5 + 1) * unsafe_fbinomial(dj5, pm456) *
            unsafe_fbinomial(j789 + 1, dj9 + 1) * unsafe_fbinomial(dj9, pm798)
    P0_de = unsafe_fbinomial(j147 + 1, dj1 + 1) * unsafe_fbinomial(dj1, div(dj1 + dj4 - dj7, 2)) *
            unsafe_fbinomial(j258 + 1, dj5 + 1) * unsafe_fbinomial(dj5, div(dj2 + dj5 - dj8, 2)) *
            unsafe_fbinomial(j369 + 1, dj9 + 1) * unsafe_fbinomial(dj9, div(dj3 + dj9 - dj6, 2))
    P0 = P0_nu / P0_de
    dtl::Int64 = max(abs(dj2 - dj6), abs(dj4 - dj8), abs(dj1 - dj9))
    dth::Int64 = min(dj2 + dj6, dj4 + dj8, dj1 + dj9)
    PABC = zero(Float64)
    for dt::Int64 = dtl:2:dth
        j19t::Int64 = div(dj1 + dj9 + dt, 2)
        j26t::Int64 = div(dj2 + dj6 + dt, 2)
        j48t::Int64 = div(dj4 + dj8 + dt, 2)
        Pt_de = unsafe_fbinomial(j19t + 1, dt + 1) * unsafe_fbinomial(dt, div(dj1 + dt - dj9, 2)) *
                unsafe_fbinomial(j26t + 1, dt + 1) * unsafe_fbinomial(dt, div(dj2 + dt - dj6, 2)) *
                unsafe_fbinomial(j48t + 1, dt + 1) * unsafe_fbinomial(dt, div(dj4 + dt - dj8, 2))
        Pt_de *= (dt + 1)^2
        xl::Int64 = max(j123, j369, j26t, j19t)
        xh::Int64 = min(pm123 + j369, pm132 + j26t, pm231 + j19t)
        At = zero(Float64)
        for x = xl:xh
            At = -At + unsafe_fbinomial(x + 1, j123 + 1) * unsafe_fbinomial(pm123, x - j369) *
                       unsafe_fbinomial(pm132, x - j26t) * unsafe_fbinomial(pm231, x - j19t)
        end
        yl::Int64 = max(j456, j26t, j258, j48t)
        yh::Int64 = min(pm456 + j26t, pm465 + j258, pm564 + j48t)
        Bt = zero(Float64)
        for y = yl:yh
            Bt = -Bt + unsafe_fbinomial(y + 1, j456 + 1) * unsafe_fbinomial(pm456, y - j26t) *
                       unsafe_fbinomial(pm465, y - j258) * unsafe_fbinomial(pm564, y - j48t)
        end
        zl::Int64 = max(j789, j19t, j48t, j147)
        zh::Int64 = min(pm789 + j19t, pm798 + j48t, pm897 + j147)
        Ct = zero(Float64)
        for z = zl:zh
            Ct = -Ct + unsafe_fbinomial(z + 1, j789 + 1) * unsafe_fbinomial(pm789, z - j19t) *
                       unsafe_fbinomial(pm798, z - j48t) * unsafe_fbinomial(pm897, z - j147)
        end
        PABC += (iphase(xh + yh + zh) * At * Bt * Ct) / Pt_de
    end
    return iphase(dth) * sqrt(P0) * PABC
end

function _fMoshinsky(N::Int, L::Int, n::Int, l::Int, n1::Int, l1::Int, n2::Int, l2::Int, Λ::Int, tanβ::Float64)
    f1 = 2 * n1 + l1
    f2 = 2 * n2 + l2
    F = 2 * N + L
    f = 2 * n + l
    f1 + f2 == F + f || return zero(Float64)

    secβ = √(1 + tanβ * tanβ)
    cosβ = 1 / secβ
    sinβ = tanβ / secβ

    nl1 = n1 + l1
    nl2 = n2 + l2
    NL = N + L
    nl = n + l
    r1 = unsafe_fbinomial(2nl1 + 1, nl1) / (unsafe_fbinomial(f1 + 2, n1) * ((nl1 + 2) << l1))
    r2 = unsafe_fbinomial(2nl2 + 1, nl2) / (unsafe_fbinomial(f2 + 2, n2) * ((nl2 + 2) << l2))
    R = unsafe_fbinomial(2NL + 1, NL) / (unsafe_fbinomial(F + 2, N) * ((NL + 2) << L))
    r = unsafe_fbinomial(2nl + 1, nl) / (unsafe_fbinomial(f + 2, n) * ((nl + 2) << l))
    pre_sum = √(r1 * r2 * R * r)
    sum = zero(Float64)
    for fa = 0:min(f1, F)
        fb = f1 - fa
        fc = F - fa
        fd = f2 - fc
        fd >= 0 || continue
        t = sinβ^(fa + fd) * cosβ^(fb + fc)
        t = t * √(unsafe_fbinomial(f1 + 2, fa + 1) * unsafe_fbinomial(f2 + 2, fc + 1) * unsafe_fbinomial(F + 2, fa + 1) * unsafe_fbinomial(f + 2, fb + 1))
        for la = (fa&0x01):2:fa
            na = div(fa - la, 2)
            nla = na + la
            ta = (((2 * la + 1) << la) * unsafe_fbinomial(fa + 1, na)) / unsafe_fbinomial(2nla + 1, nla)
            for lb = abs(l1 - la):2:min(la + l1, fb)
                nb = div(fb - lb, 2)
                nlb = nb + lb
                tb = (((2 * lb + 1) << lb) * unsafe_fbinomial(fb + 1, nb)) / unsafe_fbinomial(2nlb + 1, nlb)
                g1 = div(la + lb + l1, 2)
                CGab = unsafe_fbinomial(g1, l1) * unsafe_fbinomial(l1, g1 - la) / √(unsafe_fbinomial(2g1 + 1, 2(g1 - l1)) * unsafe_fbinomial(2l1, 2(g1 - la)))
                for lc = abs(L - la):2:min(la + L, fc)
                    nc = div(fc - lc, 2)
                    nlc = nc + lc
                    tc = (((2 * lc + 1) << lc) * unsafe_fbinomial(fc + 1, nc)) / unsafe_fbinomial(2nlc + 1, nlc)
                    G = div(la + lc + L, 2)
                    CGac = unsafe_fbinomial(G, L) * unsafe_fbinomial(L, G - la) / √(unsafe_fbinomial(2G + 1, 2(G - L)) * unsafe_fbinomial(2L, 2(G - la)))
                    ld_min = max(abs(l2 - lc), abs(l - lb))
                    ld_max = min(fd, lb + l, lc + l2)
                    for ld = ld_min:2:ld_max
                        nd = div(fd - ld, 2)
                        nld = nd + ld
                        td = (((2 * ld + 1) << ld) * unsafe_fbinomial(fd + 1, nd)) / unsafe_fbinomial(2nld + 1, nld)
                        g2 = div(lc + ld + l2, 2)
                        CGcd = unsafe_fbinomial(g2, l2) * unsafe_fbinomial(l2, g2 - lc) / √(unsafe_fbinomial(2g2 + 1, 2(g2 - l2)) * unsafe_fbinomial(2l2, 2(g2 - lc)))
                        g = div(lb + ld + l, 2)
                        CGbd = unsafe_fbinomial(g, l) * unsafe_fbinomial(l, g - lb) / √(unsafe_fbinomial(2g + 1, 2(g - l)) * unsafe_fbinomial(2l, 2(g - lb)))
                        phase = iphase(ld)
                        ninej = f9j(2 * la, 2 * lb, 2 * l1, 2 * lc, 2 * ld, 2 * l2, 2 * L, 2 * l, 2 * Λ)
                        sum = sum + phase * t * ta * tb * tc * td * CGab * CGac * CGbd * CGcd * ninej
                    end
                end
            end
        end
    end
    return pre_sum * sum
end


function _dfunc(dj::Int, dm1::Int, dm2::Int, β::Float64)
    check_jm(dj, dm1) && check_jm(dj, dm2) || return zero(Float64)
    jm1 = div(dj - dm1, 2)
    jp1 = div(dj + dm1, 2)
    jm2 = div(dj - dm2, 2)
    mm = div(dm1 + dm2, 2)
    s, c = sincos(β)
    kmin = max(0, -mm)
    kmax = min(jm1, jm2)
    sum = 0.0
    for k = kmin:kmax
        sum = -sum + unsafe_fbinomial(jm1, k) * unsafe_fbinomial(jp1, mm+k) * c^(mm + 2k) * s^(jm1 + jm2 - 2k)
    end
    sum = iphase(jm2 + kmax) * sum
    sum = sum * √(unsafe_fbinomial(dj, jm1) / unsafe_fbinomial(dj, jm2))
    return sum
end

"""
    fCG(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer)
float64 and fast CG coefficient.
"""
@inline function fCG(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer)
    return _fCG(Int64.((dj1, dj2, dj3, dm1, dm2, dm3))...)
end

"""
    fCG0(dj1::Integer, dj2::Integer, dj3::Integer)
"""
@inline function fCG0(dj1::Integer, dj2::Integer, dj3::Integer)
    return _fCG0(Int64.((dj1, dj2, dj3))...)
end

"""
    f3j(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer)
float64 and fast Wigner 3j symbol.
"""
@inline function f3j(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer)
    return iphase(dj1 + div(dj3 + dm3, 2)) * fCG(dj1, dj2, dj3, -dm1, -dm2, dm3) / sqrt(dj3 + 1)
end

"""
    f6j(dj1::Integer, dj2::Integer, dj3::Integer, dj4::Integer, dj5::Integer, dj6::Integer)
float64 and fast Wigner 6j symbol.
"""
@inline function f6j(dj1::Integer, dj2::Integer, dj3::Integer, dj4::Integer, dj5::Integer, dj6::Integer)
    return _f6j(Int64.((dj1, dj2, dj3, dj4, dj5, dj6))...)
end

"""
    Racah(dj1::Integer, dj2::Integer, dj3::Integer, dj4::Integer, dj5::Integer, dj6::Integer)
float64 and fast Racah coefficient.
"""
@inline function fRacah(dj1::Integer, dj2::Integer, dj3::Integer, dj4::Integer, dj5::Integer, dj6::Integer)
    return iphase(div(dj1 + dj2 + dj3 + dj4, 2)) * f6j(dj1, dj2, dj5, dj4, dj3, dj6)
end

"""
    f9j(dj1::Integer, dj2::Integer, dj3::Integer,
        dj4::Integer, dj5::Integer, dj6::Integer,
        dj7::Integer, dj8::Integer, dj9::Integer)
float64 and fast Wigner 9j symbol.
"""
@inline function f9j(dj1::Integer, dj2::Integer, dj3::Integer,
    dj4::Integer, dj5::Integer, dj6::Integer,
    dj7::Integer, dj8::Integer, dj9::Integer)
    return _f9j(Int64.((dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9))...)
end

"""
    fMoshinsky(N::Integer, L::Integer, n::Integer, l::Integer, n1::Integer, l1::Integer, n2::Integer, l2::Integer, tanβ::Float64 = 1.0)
float64 and fast Moshinsky bracket.
"""
@inline function fMoshinsky(N::Integer, L::Integer, n::Integer, l::Integer, n1::Integer, l1::Integer, n2::Integer, l2::Integer, Λ::Integer, tanβ::Float64=1.0)
    return _fMoshinsky(Int.((N, L, n, l, n1, l1, n2, l2, Λ))..., tanβ)
end

"""
    dfunc(dj::Integer, dm1::Integer, dm2::Integer, β::Float64)
Wigner-d function.
"""
@inline function dfunc(dj::Integer, dm1::Integer, dm2::Integer, β::Float64)
    return _dfunc(Int.((dj, dm1, dm2))..., β)
end