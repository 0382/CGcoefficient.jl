_fbinomial_data = Float64[1.0]
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

function _extend_fbinomial_data(n::Int)
    if n > get_fbinomial_nmax()
        old_data = copy(get_fbinomial_data())
        resize!(get_fbinomial_data(), binomial_data_size(n))
        copyto!(get_fbinomial_data(), old_data)
        t = BigInt()
        for m = get_fbinomial_nmax()+1:n
            for k = 0:div(m, 2)
                get_fbinomial_data()[binomial_index(m, k)] = _bigbin(t, m, k)
            end
            global _fbinomial_nmax = get_fbinomial_nmax() + 1
        end
        global _fbinomial_nmax = n
    end
    nothing
end

"""
    wigner_init_float(n::Integer, mode::AbstractString, rank::Integer)
We calculate the Wigner Symbols with binomials. We will store some binomials first, then when we need one binomial, we just load it. In this package, the functions starts with `f` character dependent on the stored binomials. The `wigner_init_float` function is used to extend the range of the stored binomial table.

The parameters means:
- `mode`:
    - `"Jmax"`: set the maximum angular momentum in the system.
    - `"2bjmax"`: set the maximum two-body coupled angular momentum in the system.
    - `"nmax"`: directly set the maximum binomial number.
    - `"Moshinsky"`: set the maximum `N = 2n+l` in the Moshinsky bracket.
- `n`: the maximum angular momentum or the maximum binomial number, dependent on the `mode`.
- `rank`:
    - `3`: calculate only CG coefficients and 3j symbols.
    - `6`: calculate CG coefficients, 3j symbols and 6j symbols.
    - `9`: calculate CG coefficients, 3j symbols, 6j symbols and 9j symbols.

`"Jmax"` means the global maximum angular momentum, for every parameters. It is always safe with out any assumption.

The `"2bjmax"` mode means your calculation only consider two-body coupling, and no three-body coupling. This mode assumes that in all these coefficients, at least one of the angular momentun is just a single particle angular momentum. With this assumption, `"2bjmax"` mode will use less memory than `"Jmax"` mode.

The `"Moshinsky"` mode means you want to calculate the Moshinsky brackets. The `n` parameter is the maximum `N = 2n+l` in the Moshinsky bracket. The `rank` parameter is ignored.

The `"nmax"` mode directly set `nmax`, and the `rank` parameter is ignored. 

For example
```julia
wigner_init_float(21, "Jmax", 6)
```
means you calculate CG and 6j symbols, and donot calculate 9j symbol. The maximum angular momentum in your system is 21.

The exact values of the maximum `nmax` for different `rank` are shown in the table below ([details](https://0382.github.io/CGcoefficient.jl/stable/wigner/#Estimate-the-capacity)).

|                                       |    Calculate range    |   CG & 3j   | 6j & Racah  |     9j      |
| :-----------------------------------: | :-------------------: | :---------: | :---------: | :---------: |
|          meaning of `type`            | `type`\\\\`rank`      |      3      |      6      |      9      |
|        max angular momentum           |    `"Jmax"`           |  `3Jmax+1`  |  `4Jmax+1`  |  `5Jmax+1`  |
| max two-body coupled angular momentum |   `"2bjmax"`          |  `2jmax+1`  |  `3jmax+1`  |  `4jmax+1`  |
| max ``N = 2n+l`` for Moshinsky bracket|    `"Moshinsky"`      |  `6Nmax+1`  |  `6Nmax+1`  |  `6Nmax+1`  |
|            max binomial               |    `"nmax"`           |   `nmax`    |   `namx`    |   `nmax`    |

You do not need to rememmber those values in the table. You just need to find the maximum angular momentum in you canculation, then call the function.

The `wigner_init_float` function is **not** thread safe, so you should call it before you start your calculation.
"""
function wigner_init_float(n::Integer, mode::AbstractString, rank::Integer)
    if mode == "Jmax"
        rank == 3 && _extend_fbinomial_data(3 * n + 1)
        rank == 6 && _extend_fbinomial_data(4 * n + 1)
        rank == 9 && _extend_fbinomial_data(5 * n + 1)
        if rank != 3 && rank != 6 && rank != 9
            throw(ArgumentError("invalid rank $rank"))
        end
    elseif mode == "2bjmax"
        rank == 3 && _extend_fbinomial_data(2 * n + 1)
        rank == 6 && _extend_fbinomial_data(3 * n + 1)
        rank == 9 && _extend_fbinomial_data(4 * n + 1)
        if rank != 3 && rank != 6 && rank != 9
            throw(ArgumentError("invalid rank $rank"))
        end
    elseif mode == "nmax"
        _extend_fbinomial_data(n)
    elseif mode == "Moshinsky"
        _extend_fbinomial_data(6 * n + 1)
    else
        throw(ArgumentError("invalid mode $mode"))
    end
    nothing
end


# basic CG coefficient calculation function with Float64 return value
# unlike the SqrtRational version, this function gives out zero when the parameters are not valid
function _fCG(dj1::Int, dj2::Int, dj3::Int, dm1::Int, dm2::Int, dm3::Int)
    check_CG(dj1, dj2, dj3, dm1, dm2, dm3) || return zero(Float64)
    J::Int = div(dj1 + dj2 + dj3, 2)
    Jm1::Int = J - dj1
    Jm2::Int = J - dj2
    Jm3::Int = J - dj3
    j1mm1::Int = div(dj1 - dm1, 2)
    j2mm2::Int = div(dj2 - dm2, 2)
    j3mm3::Int = div(dj3 - dm3, 2)
    j2pm2::Int = div(dj2 + dm2, 2)
    A = sqrt(unsafe_fbinomial(dj1, Jm2) * unsafe_fbinomial(dj2, Jm3) / (
        unsafe_fbinomial(J + 1, Jm3) * unsafe_fbinomial(dj1, j1mm1) *
        unsafe_fbinomial(dj2, j2mm2) * unsafe_fbinomial(dj3, j3mm3)
    ))
    B = zero(Float64)
    low::Int = max(zero(Int), j1mm1 - Jm2, j2pm2 - Jm1)
    high::Int = min(Jm3, j1mm1, j2pm2)
    for z in low:high
        B = -B + unsafe_fbinomial(Jm3, z) * unsafe_fbinomial(Jm2, j1mm1 - z) * unsafe_fbinomial(Jm1, j2pm2 - z)
    end
    return iphase(high) * A * B
end

function _f3j(dj1::Int, dj2::Int, dj3::Int, dm1::Int, dm2::Int, dm3::Int)
    check_3j(dj1, dj2, dj3, dm1, dm2, dm3) || return zero(Float64)
    J::Int = div(dj1 + dj2 + dj3, 2)
    Jm1::Int = J - dj1
    Jm2::Int = J - dj2
    Jm3::Int = J - dj3
    j1mm1::Int = div(dj1 - dm1, 2)
    j2mm2::Int = div(dj2 - dm2, 2)
    j3mm3::Int = div(dj3 - dm3, 2)
    j1pm1::Int = div(dj1 + dm1, 2)
    A = sqrt(unsafe_fbinomial(dj1, Jm2) * unsafe_fbinomial(dj2, Jm1) / (
        (J+1) * unsafe_fbinomial(J, Jm3) * unsafe_fbinomial(dj1, j1mm1) *
        unsafe_fbinomial(dj2, j2mm2) * unsafe_fbinomial(dj3, j3mm3)
    ))
    B = zero(Float64)
    low::Int = max(zero(Int), j1pm1 - Jm2, j2mm2 - Jm1)
    high::Int = min(Jm3, j1pm1, j2mm2)
    for z in low:high
        B = -B + unsafe_fbinomial(Jm3, z) * unsafe_fbinomial(Jm2, j1pm1 - z) * unsafe_fbinomial(Jm1, j2mm2 - z)
    end
    return iphase(dj1 + div(dj3 + dm3, 2) + high) * A * B
end

function _fCG0(j1::Int, j2::Int, j3::Int)
    check_couple(2j1, 2j2, 2j3) || return zero(Float64)
    J = j1 + j2 + j3
    isodd(J) && return zero(Float64)
    g = div(J, 2)
    return iphase(g - j3) * unsafe_fbinomial(g, j3) * unsafe_fbinomial(j3, g - j1) / sqrt(unsafe_fbinomial(J + 1, 2j3 + 1) * unsafe_fbinomial(2j3, J - 2j1))
end

function _f6j(dj1::Int, dj2::Int, dj3::Int, dj4::Int, dj5::Int, dj6::Int)
    check_6j(dj1, dj2, dj3, dj4, dj5, dj6) || return zero(Float64)
    j123::Int = div(dj1 + dj2 + dj3, 2)
    j156::Int = div(dj1 + dj5 + dj6, 2)
    j426::Int = div(dj4 + dj2 + dj6, 2)
    j453::Int = div(dj4 + dj5 + dj3, 2)
    jpm123::Int = div(dj1 + dj2 - dj3, 2)
    jpm132::Int = div(dj1 + dj3 - dj2, 2)
    jpm231::Int = div(dj2 + dj3 - dj1, 2)
    jpm156::Int = div(dj1 + dj5 - dj6, 2)
    jpm426::Int = div(dj4 + dj2 - dj6, 2)
    jpm453::Int = div(dj4 + dj5 - dj3, 2)
    A = sqrt(unsafe_fbinomial(j123 + 1, dj1 + 1) * unsafe_fbinomial(dj1, jpm123) / (
        unsafe_fbinomial(j156 + 1, dj1 + 1) * unsafe_fbinomial(dj1, jpm156) *
        unsafe_fbinomial(j453 + 1, dj4 + 1) * unsafe_fbinomial(dj4, jpm453) *
        unsafe_fbinomial(j426 + 1, dj4 + 1) * unsafe_fbinomial(dj4, jpm426)
    ))
    B = zero(Float64)
    low::Int = max(j123, j453, j426, j156)
    high::Int = min(jpm123 + j453, jpm132 + j426, jpm231 + j156)
    for x = low:high
        B = -B + unsafe_fbinomial(x + 1, j123 + 1) * unsafe_fbinomial(jpm123, x - j453) *
                 unsafe_fbinomial(jpm132, x - j426) * unsafe_fbinomial(jpm231, x - j156)
    end
    return iphase(high) * A * B / (dj4 + 1)
end

function _f9j(dj1::Int, dj2::Int, dj3::Int,
    dj4::Int, dj5::Int, dj6::Int,
    dj7::Int, dj8::Int, dj9::Int)
    check_9j(dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9) || return zero(Float64)
    j123::Int = div(dj1 + dj2 + dj3, 2)
    j456::Int = div(dj4 + dj5 + dj6, 2)
    j789::Int = div(dj7 + dj8 + dj9, 2)
    j147::Int = div(dj1 + dj4 + dj7, 2)
    j258::Int = div(dj2 + dj5 + dj8, 2)
    j369::Int = div(dj3 + dj6 + dj9, 2)
    pm123::Int = div(dj1 + dj2 - dj3, 2)
    pm132::Int = div(dj1 + dj3 - dj2, 2)
    pm231::Int = div(dj2 + dj3 - dj1, 2)
    pm456::Int = div(dj4 + dj5 - dj6, 2)
    pm465::Int = div(dj4 + dj6 - dj5, 2)
    pm564::Int = div(dj5 + dj6 - dj4, 2)
    pm789::Int = div(dj7 + dj8 - dj9, 2)
    pm798::Int = div(dj7 + dj9 - dj8, 2)
    pm897::Int = div(dj8 + dj9 - dj7, 2)
    # value in sqrt
    P0_nu = unsafe_fbinomial(j123 + 1, dj1 + 1) * unsafe_fbinomial(dj1, pm123) *
            unsafe_fbinomial(j456 + 1, dj5 + 1) * unsafe_fbinomial(dj5, pm456) *
            unsafe_fbinomial(j789 + 1, dj9 + 1) * unsafe_fbinomial(dj9, pm798)
    P0_de = unsafe_fbinomial(j147 + 1, dj1 + 1) * unsafe_fbinomial(dj1, div(dj1 + dj4 - dj7, 2)) *
            unsafe_fbinomial(j258 + 1, dj5 + 1) * unsafe_fbinomial(dj5, div(dj2 + dj5 - dj8, 2)) *
            unsafe_fbinomial(j369 + 1, dj9 + 1) * unsafe_fbinomial(dj9, div(dj3 + dj9 - dj6, 2))
    P0 = P0_nu / P0_de
    dtl::Int = max(abs(dj2 - dj6), abs(dj4 - dj8), abs(dj1 - dj9))
    dth::Int = min(dj2 + dj6, dj4 + dj8, dj1 + dj9)
    PABC = zero(Float64)
    for dt::Int = dtl:2:dth
        j19t::Int = div(dj1 + dj9 + dt, 2)
        j26t::Int = div(dj2 + dj6 + dt, 2)
        j48t::Int = div(dj4 + dj8 + dt, 2)
        Pt_de = unsafe_fbinomial(j19t + 1, dt + 1) * unsafe_fbinomial(dt, div(dj1 + dt - dj9, 2)) *
                unsafe_fbinomial(j26t + 1, dt + 1) * unsafe_fbinomial(dt, div(dj2 + dt - dj6, 2)) *
                unsafe_fbinomial(j48t + 1, dt + 1) * unsafe_fbinomial(dt, div(dj4 + dt - dj8, 2))
        Pt_de *= (dt + 1)^2
        xl::Int = max(j123, j369, j26t, j19t)
        xh::Int = min(pm123 + j369, pm132 + j26t, pm231 + j19t)
        At = zero(Float64)
        for x = xl:xh
            At = -At + unsafe_fbinomial(x + 1, j123 + 1) * unsafe_fbinomial(pm123, x - j369) *
                       unsafe_fbinomial(pm132, x - j26t) * unsafe_fbinomial(pm231, x - j19t)
        end
        yl::Int = max(j456, j26t, j258, j48t)
        yh::Int = min(pm456 + j26t, pm465 + j258, pm564 + j48t)
        Bt = zero(Float64)
        for y = yl:yh
            Bt = -Bt + unsafe_fbinomial(y + 1, j456 + 1) * unsafe_fbinomial(pm456, y - j26t) *
                       unsafe_fbinomial(pm465, y - j258) * unsafe_fbinomial(pm564, y - j48t)
        end
        zl::Int = max(j789, j19t, j48t, j147)
        zh::Int = min(pm789 + j19t, pm798 + j48t, pm897 + j147)
        Ct = zero(Float64)
        for z = zl:zh
            Ct = -Ct + unsafe_fbinomial(z + 1, j789 + 1) * unsafe_fbinomial(pm789, z - j19t) *
                       unsafe_fbinomial(pm798, z - j48t) * unsafe_fbinomial(pm897, z - j147)
        end
        PABC += (iphase(xh + yh + zh) * At * Bt * Ct) / Pt_de
    end
    return iphase(dth) * sqrt(P0) * PABC
end

function _fm9j(j1::Int, j2::Int, j3::Int, j4::Int, j5::Int, j6::Int, j7::Int, j8::Int, j9::Int)::Float64
    j123::Int = j1 + j2 + j3
    j456::Int = j4 + j5 + j6
    j789::Int = j7 + j8 + j9
    j147::Int = j1 + j4 + j7
    j258::Int = j2 + j5 + j8
    j369::Int = j3 + j6 + j9
    pm123::Int = j1 + j2 - j3
    pm132::Int = j1 + j3 - j2
    pm231::Int = j2 + j3 - j1
    pm456::Int = j4 + j5 - j6
    pm465::Int = j4 + j6 - j5
    pm564::Int = j5 + j6 - j4
    pm789::Int = j7 + j8 - j9
    pm798::Int = j7 + j9 - j8
    pm897::Int = j8 + j9 - j7
    sum = zero(Float64)
    tl::Int = max(abs(j2 - j6), abs(j4 - j8), abs(j1 - j9))
    th::Int = min(j2 + j6, j4 + j8, j1 + j9)
    for t::Int = tl:th
        j19t::Int = j1 + j9 + t
        j26t::Int = j2 + j6 + t
        j48t::Int = j4 + j8 + t
        dt::Int = 2t
        Pt_de = unsafe_fbinomial(j19t + 1, dt + 1) * unsafe_fbinomial(dt, j1 + t - j9) *
                unsafe_fbinomial(j26t + 1, dt + 1) * unsafe_fbinomial(dt, j2 + t - j6) *
                unsafe_fbinomial(j48t + 1, dt + 1) * unsafe_fbinomial(dt, j4 + t - j8)
        Pt_de *= (dt + 1)^2
        xl::Int = max(j123, j369, j26t, j19t)
        xh::Int = min(pm123 + j369, pm132 + j26t, pm231 + j19t)
        At = zero(Float64)
        for x::Int = xl:xh
            At = -At + unsafe_fbinomial(x + 1, j123 + 1) * unsafe_fbinomial(pm123, x - j369) *
                       unsafe_fbinomial(pm132, x - j26t) * unsafe_fbinomial(pm231, x - j19t)
        end
        yl::Int = max(j456, j26t, j258, j48t)
        yh::Int = min(pm456 + j26t, pm465 + j258, pm564 + j48t)
        Bt = zero(Float64)
        for y::Int = yl:yh
            Bt = -Bt + unsafe_fbinomial(y + 1, j456 + 1) * unsafe_fbinomial(pm456, y - j26t) *
                       unsafe_fbinomial(pm465, y - j258) * unsafe_fbinomial(pm564, y - j48t)
        end
        zl::Int = max(j789, j19t, j48t, j147)
        zh::Int = min(pm789 + j19t, pm798 + j48t, pm897 + j147)
        Ct = zero(Float64)
        for z::Int = zl:zh
            Ct = -Ct + unsafe_fbinomial(z + 1, j789 + 1) * unsafe_fbinomial(pm789, z - j19t) *
                       unsafe_fbinomial(pm798, z - j48t) * unsafe_fbinomial(pm897, z - j147)
        end
        sum += iphase(xh + yh + zh) * At * Bt * Ct / Pt_de
    end
    return sum
end

function _fMoshinsky(N::Int, L::Int, n::Int, l::Int, n1::Int, l1::Int, n2::Int, l2::Int, Λ::Int, tanβ::Float64)
    check_Moshinsky(N, L, n, l, n1, l1, n2, l2, Λ) || return zero(Float64)
    e1 = 2 * n1 + l1
    e2 = 2 * n2 + l2
    E = 2 * N + L
    e = 2 * n + l
    χ = e1 + e2
    nl1 = n1 + l1
    nl2 = n2 + l2
    NL = N + L
    nl = n + l

    cosβ = 1 / sqrt(1 + tanβ^2)
    sinβ = tanβ * cosβ
    pre = unsafe_fbinomial(χ + 2, e1 + 1) / unsafe_fbinomial(χ + 2, E + 1)
    pre *= unsafe_fbinomial(L + l + Λ + 1, 2Λ + 1) * unsafe_fbinomial(2Λ, Λ + L - l)
    pre /= unsafe_fbinomial(l1 + l2 + Λ + 1, 2Λ + 1) * unsafe_fbinomial(2Λ, Λ + l1 - l2)

    pre *= (2l1 + 1) * unsafe_fbinomial(2nl1 + 1, nl1) / (unsafe_fbinomial(e1 + 1, n1) * 2^l1)
    pre *= (2l2 + 1) * unsafe_fbinomial(2nl2 + 1, nl2) / (unsafe_fbinomial(e2 + 1, n2) * 2^l2)
    pre *= (2L + 1) * unsafe_fbinomial(2NL + 1, NL) / (unsafe_fbinomial(E + 1, N) * 2^L)
    pre *= (2l + 1) * unsafe_fbinomial(2nl + 1, nl) / (unsafe_fbinomial(e + 1, n) * 2^l)
    pre = sqrt(pre) / ((e1 + 2) * (e2 + 2))

    sum = zero(Float64)
    for ea = 0:min(e1, E)
        eb = e1 - ea
        ec = E - ea
        ed = e2 - ec
        ed < 0 && continue
        tfa = sinβ^(ea + ed) * cosβ^(eb + ec) * unsafe_fbinomial(e1 + 2, ea + 1) * unsafe_fbinomial(e2 + 2, ec + 1)
        for la = rem(ea, 2):2:ea
            na = div(ea - la, 2)
            nla = na + la
            ta = 2^la * (2la + 1) * unsafe_fbinomial(ea + 1, na) / unsafe_fbinomial(2nla + 1, nla)
            ta = tfa * ta
            for lb = abs(l1 - la):2:min(l1 + la, eb)
                nb = div(eb - lb, 2)
                nlb = nb + lb
                tb = 2^lb * (2lb + 1) * unsafe_fbinomial(eb + 1, nb) / unsafe_fbinomial(2nlb + 1, nlb)
                g1 = div(la + lb + l1, 2)
                t1 = unsafe_fbinomial(g1, l1) * unsafe_fbinomial(l1, g1 - la)
                tb = ta * tb * t1
                for lc = abs(L - la):2:min(L + la, ec)
                    nc = div(ec - lc, 2)
                    nlc = nc + lc
                    tc = 2^lc * (2lc + 1) * unsafe_fbinomial(ec + 1, nc) / unsafe_fbinomial(2nlc + 1, nlc)
                    g3 = div(la + lc + L, 2)
                    t3 = unsafe_fbinomial(g3, L) * unsafe_fbinomial(L, g3 - la)
                    Δ3 = (2L + 1) * unsafe_fbinomial(la + lc + L + 1, 2L + 1) * unsafe_fbinomial(2L, L + la - lc)
                    tc = tb * tc * t3 / Δ3
                    ldmin = max(abs(l2 - lc), abs(l - lb))
                    ldmax = min(ed, l2 + lc, l + lb)
                    for ld = ldmin:2:ldmax
                        nd = div(ed - ld, 2)
                        nld = nd + ld
                        td = 2^ld * (2ld + 1) * unsafe_fbinomial(ed + 1, nd) / unsafe_fbinomial(2nld + 1, nld)
                        g2 = div(lc + ld + l2, 2)
                        t2 = unsafe_fbinomial(g2, l2) * unsafe_fbinomial(l2, g2 - lc)
                        g4 = div(lb + ld + l, 2)
                        t4 = unsafe_fbinomial(g4, l) * unsafe_fbinomial(l, g4 - lb)
                        Δ4 = (2l + 1) * unsafe_fbinomial(lb + ld + l + 1, 2l + 1) * unsafe_fbinomial(2l, l + lb - ld)
                        td = tc * td * t2 * t4 / Δ4
                        m9j = _fm9j(la, lb, l1, lc, ld, l2, L, l, Λ)
                        sum += iphase(ld) * td * m9j
                    end
                end
            end
        end
    end
    return pre * sum
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
        sum = -sum + unsafe_fbinomial(jm1, k) * unsafe_fbinomial(jp1, mm + k) * c^(mm + 2k) * s^(jm1 + jm2 - 2k)
    end
    sum = iphase(jm2 + kmax) * sum
    sum = sum * √(unsafe_fbinomial(dj, jm1) / unsafe_fbinomial(dj, jm2))
    return sum
end

function _flsjj(l1::Int, l2::Int, dj1::Int, dj2::Int, L::Int, S::Int, J::Int)::Float64
    check_couple(2l1, 2l2, 2L) || return zero(Float64)
    check_couple(dj1, dj2, 2J) || return zero(Float64)
    LSJcase = 0
    if S == 0
        LSJcase = 1
    elseif S == 1
        (L == J - 1) && (LSJcase = 2)
        (J == L) && (LSJcase = 3)
        (L == J + 1) && (LSJcase = 4)
    end
    iszero(LSJcase) && return zero(Float64)
    m = l1 - l2
    p = l1 + l2
    δ1 = dj1 - 2l1
    δ2 = dj2 - 2l2
    d = 2*(2l1+1)*(2l2+1)
    r = zero(Float64)
    if δ1 == 1 && δ2 == 1
        if LSJcase == 1
            r = (J + p + 2) * (p + 1 - J) / d
        elseif LSJcase == 2
            r = (J + m) * (J - m) / (J * (2J + 1))
            r *= (J + p + 1) * (J + p + 2) / d
        elseif LSJcase == 3
            m == 0 && return zero(Float64)
            r = m * abs(m) / (J * (J + 1))
            r *= (J + p + 2) * (p + 1 - J) / d
        elseif LSJcase == 4
            r = -(L + m) * (L - m) / (L * (2J + 1))
            r *= (p - J) * (p + 1 - J) / d
        end
    elseif δ1 == 1 && δ2 == -1
        if LSJcase == 1
            r = (J + m + 1) * (J - m) / d
        elseif LSJcase == 2
            r = -(p + 1 + J) * (p + 1 - J) / (J * (2J + 1))
            r *= (J + m + 1) * (J + m) / d
        elseif LSJcase == 3
            r = (p + 1) * (p + 1) / (J * (J + 1))
            r *= (J + m + 1) * (J - m) / d
        elseif LSJcase == 4
            r = -(p + 1 + L) * (p + 1 - L) / (L * (2J + 1))
            r *= (J - m + 1) * (J - m) / d
        end
    elseif δ1 == -1 && δ2 == 1
        if LSJcase == 1
            r = -(J + m) * (J - m + 1) / d
        elseif LSJcase == 2
            r = (p + 1 + J) * (p + 1 - J) / (J * (2J + 1))
            r *= (J - m) * (J - m + 1) / d
        elseif LSJcase == 3
            r = (p + 1) * (p + 1) / (J * (J + 1))
            r *= (J + m) * (J - m + 1) / d
        elseif LSJcase == 4
            r = (p + 1 + L) * (p + 1 - L) / (L * (2J + 1))
            r *= (J + m + 1) * (J + m) / d
        end
    elseif δ1 == -1 && δ2 == -1
        if LSJcase == 1
            r = (J + p + 1) * (p - J) / d
        elseif LSJcase == 2
            r = -(J + m) * (J - m) / (J * (2J + 1))
            r *= (p - J) * (p + 1 - J) / d
        elseif LSJcase == 3
            m == 0 && return zero(Float64)
            r = -m * abs(m) / (J * (J + 1))
            r *= (J + p + 1) * (p - J) / d
        elseif LSJcase == 4
            r = (L + m) * (L - m) / (L * (2J + 1))
            r *= (J + p + 1) * (J + p + 2) / d
        end
    end
    return copysign(sqrt(abs(r)), r)
end

"""
    fCG(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer)
float64 and fast CG coefficient.
"""
@inline function fCG(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer)
    return _fCG(Int.((dj1, dj2, dj3, dm1, dm2, dm3))...)
end

"""
    fCG0(dj1::Integer, dj2::Integer, dj3::Integer)
float64 and fast CG coefficient for `m1 == m2 == m3 == 0`.
"""
@inline function fCG0(dj1::Integer, dj2::Integer, dj3::Integer)
    return _fCG0(Int.((dj1, dj2, dj3))...)
end

"""
    fCGspin(ds1::Integer, ds2::Integer, S::Integer)
float64 and fast CG coefficient for two spin-1/2 coupling.
"""
@inline function fCGspin(ds1::Integer, ds2::Integer, S::Integer)
    unsigned(S) > 1 && return zero(Float64)
    (abs(ds1) != 1 || abs(ds2) != 1) && return zero(Float64)
    if iszero(S)
        return ds1 == ds2 ? 0.0 : copysign(1 / √2, ds1)
    else # S = 1
        return ds1 == ds2 ? 1.0 : 1 / √2
    end
end

# i1 = ds1 > 0, i2 = ds2 > 0, i3 = S > 0, is = (S12 + dS) / 2
# (i1, i2, i3, is) -> 3 * (4 * i1 + 2 * i2 + i3) + is
const CG3spin_data = [
    0.0,                 # (-1/2, -1/2, -1/2,   0,  1/2) -> 0
    0.0,                 # (-1/2, -1/2, -1/2,   1,  1/2) -> 0
    1.0,                 # (-1/2, -1/2, -1/2,   1,  3/2) -> 1
    0.0,                 # (-1/2, -1/2,  1/2,   0,  1/2) -> 0
    -0.816496580927726,  # (-1/2, -1/2,  1/2,   1,  1/2) -> -sqrt(2/3)
    0.5773502691896257,  # (-1/2, -1/2,  1/2,   1,  3/2) -> sqrt(1/3)
    -0.7071067811865476, # (-1/2,  1/2, -1/2,   0,  1/2) -> -sqrt(1/2)
    0.408248290463863,   # (-1/2,  1/2, -1/2,   1,  1/2) -> sqrt(1/6)
    0.5773502691896257,  # (-1/2,  1/2, -1/2,   1,  3/2) -> sqrt(1/3)
    -0.7071067811865476, # (-1/2,  1/2,  1/2,   0,  1/2) -> -sqrt(1/2)
    -0.408248290463863,  # (-1/2,  1/2,  1/2,   1,  1/2) -> -sqrt(1/6)
    0.5773502691896257,  # (-1/2,  1/2,  1/2,   1,  3/2) -> sqrt(1/3)
    0.7071067811865476,  # ( 1/2, -1/2, -1/2,   0,  1/2) -> sqrt(1/2)
    0.408248290463863,   # ( 1/2, -1/2, -1/2,   1,  1/2) -> sqrt(1/6)
    0.5773502691896257,  # ( 1/2, -1/2, -1/2,   1,  3/2) -> sqrt(1/3)
    0.7071067811865476,  # ( 1/2, -1/2,  1/2,   0,  1/2) -> sqrt(1/2)
    -0.408248290463863,  # ( 1/2, -1/2,  1/2,   1,  1/2) -> -sqrt(1/6)
    0.5773502691896257,  # ( 1/2, -1/2,  1/2,   1,  3/2) -> sqrt(1/3)
    0.0,                 # ( 1/2,  1/2, -1/2,   0,  1/2) -> 0
    0.816496580927726,   # ( 1/2,  1/2, -1/2,   1,  1/2) -> sqrt(2/3)
    0.5773502691896257,  # ( 1/2,  1/2, -1/2,   1,  3/2) -> sqrt(1/3)
    0.0,                 # ( 1/2,  1/2,  1/2,   0,  1/2) -> 0
    0.0,                 # ( 1/2,  1/2,  1/2,   1,  1/2) -> 0
    1.0                  # ( 1/2,  1/2,  1/2,   1,  3/2) -> 1
]

"""
    fCG3spin(ds1::Integer, ds2::Integer, ds3::Integer, S12::Integer, dS::Integer)
float64 and fast CG coefficient for three spin-1/2 coupling: <S12,M12|1/2,m1;1/2,m2><S,M|S12,M12;1/2,m3>
"""
function fCG3spin(ds1::Integer, ds2::Integer, ds3::Integer, S12::Integer, dS::Integer)::Float64
    unsigned(S12) > 1 && return zero(Float64)
    (S12 == 0 && dS != 1) && return zero(Float64)
    (S12 == 1 && (dS != 1 && dS != 3)) && return zero(Float64)
    (abs(ds1) != 1 || abs(ds2) != 1 || abs(ds3) != 1) && return zero(Float64)
    return CG3spin_data[3 * (2ds1 + ds2 + div(ds3 + 1, 2)) + div(S12 + dS, 2) + 10]
end


"""
    f3j(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer)
float64 and fast Wigner 3j symbol.
"""
@inline function f3j(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer)
    return _f3j(Int.((dj1, dj2, dj3, dm1, dm2, dm3))...)
end

"""
    f6j(dj1::Integer, dj2::Integer, dj3::Integer, dj4::Integer, dj5::Integer, dj6::Integer)
float64 and fast Wigner 6j symbol.
"""
@inline function f6j(dj1::Integer, dj2::Integer, dj3::Integer, dj4::Integer, dj5::Integer, dj6::Integer)
    return _f6j(Int.((dj1, dj2, dj3, dj4, dj5, dj6))...)
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
    return _f9j(Int.((dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9))...)
end

"""
    fnorm9j(dj1::Integer, dj2::Integer, dj3::Integer,
            dj4::Integer, dj5::Integer, dj6::Integer,
            dj7::Integer, dj8::Integer, dj9::Integer)
float64 and fast normalized Wigner 9j symbol.
"""
@inline function fnorm9j(dj1::Integer, dj2::Integer, dj3::Integer,
    dj4::Integer, dj5::Integer, dj6::Integer,
    dj7::Integer, dj8::Integer, dj9::Integer)
    return sqrt((dj3 + 1.0) * (dj6 + 1.0) * (dj7 + 1.0) * (dj8 + 1.0)) * f9j(dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9)
end

"""
    lsjj(l1::Integer, l2::Integer, dj1::Integer, dj2::Integer, L::Integer, S::Integer, J::Integer)
float64 and fast lsjj coefficient.
"""
@inline function flsjj(l1::Integer, l2::Integer, dj1::Integer, dj2::Integer, L::Integer, S::Integer, J::Integer)
    return _flsjj(Int.((l1, l2, dj1, dj2, L, S, J))...)
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
Wigner d-function.
"""
@inline function dfunc(dj::Integer, dm1::Integer, dm2::Integer, β::Float64)
    return _dfunc(Int.((dj, dm1, dm2))..., β)
end