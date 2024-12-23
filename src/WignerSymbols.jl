# This file contains the core functions of WignerSymbols and CG-coefficient.

# basic CG coefficient calculation function
function _dCG(dj1::BigInt, dj2::BigInt, dj3::BigInt, dm1::BigInt, dm2::BigInt, dm3::BigInt)
    check_CG(dj1, dj2, dj3, dm1, dm2, dm3) || return zero(SqrtRational)
    J::BigInt = div(dj1 + dj2 + dj3, 2)
    Jm1::BigInt = J - dj1
    Jm2::BigInt = J - dj2
    Jm3::BigInt = J - dj3
    j1mm1::BigInt = div(dj1 - dm1, 2)
    j2mm2::BigInt = div(dj2 - dm2, 2)
    j3mm3::BigInt = div(dj3 - dm3, 2)
    j2pm2::BigInt = div(dj2 + dm2, 2)
    # value in sqrt
    A = (binomial(dj1, Jm2) * binomial(dj2, Jm3)) // (
        binomial(J + 1, Jm3) * binomial(dj1, j1mm1) * binomial(dj2, j2mm2) * binomial(dj3, j3mm3)
    )
    B::BigInt = zero(BigInt)
    low::BigInt = max(zero(BigInt), j1mm1 - Jm2, j2pm2 - Jm1)
    high::BigInt = min(Jm3, j1mm1, j2pm2)
    for z in low:high
        B = -B + binomial(Jm3, z) * binomial(Jm2, j1mm1 - z) * binomial(Jm1, j2pm2 - z)
    end
    return SqrtRational(iphase(high) * B, A)
end

# spaecial case: m1 == m2 == m3 == 0
function _CG0(j1::BigInt, j2::BigInt, j3::BigInt)
    check_couple(2j1, 2j2, 2j3) || return zero(SqrtRational)
    J = j1 + j2 + j3
    isodd(J) && return zero(SqrtRational)
    g = div(J, 2)
    return SqrtRational(iphase(g - j3) * binomial(g, j3) * binomial(j3, g - j1), big(1) // (binomial(J + 1, 2j3 + 1) * binomial(2j3, J - 2j1)))
end

# use CG coefficient to calculate 3j-symbol
function _d3j(dj1::BigInt, dj2::BigInt, dj3::BigInt, dm1::BigInt, dm2::BigInt, dm3::BigInt)
    iphase(dj1 + div(dj3 + dm3, 2)) * exact_sqrt(1 // (dj3 + 1)) * _dCG(dj1, dj2, dj3, -dm1, -dm2, dm3)
end

# basic 6j-symbol calculation funciton
function _d6j(dj1::BigInt, dj2::BigInt, dj3::BigInt, dj4::BigInt, dj5::BigInt, dj6::BigInt)
    check_6j(dj1, dj2, dj3, dj4, dj5, dj6) || return zero(SqrtRational)
    j123::BigInt = div(dj1 + dj2 + dj3, 2)
    j156::BigInt = div(dj1 + dj5 + dj6, 2)
    j426::BigInt = div(dj4 + dj2 + dj6, 2)
    j453::BigInt = div(dj4 + dj5 + dj3, 2)
    jpm123::BigInt = div(dj1 + dj2 - dj3, 2)
    jpm132::BigInt = div(dj1 + dj3 - dj2, 2)
    jpm231::BigInt = div(dj2 + dj3 - dj1, 2)
    jpm156::BigInt = div(dj1 + dj5 - dj6, 2)
    jpm426::BigInt = div(dj4 + dj2 - dj6, 2)
    jpm453::BigInt = div(dj4 + dj5 - dj3, 2)
    # value in sqrt
    A = (binomial(j123 + 1, dj1 + 1) * binomial(dj1, jpm123)) // (
        binomial(j156 + 1, dj1 + 1) * binomial(dj1, jpm156) *
        binomial(j453 + 1, dj4 + 1) * binomial(dj4, jpm453) *
        binomial(j426 + 1, dj4 + 1) * binomial(dj4, jpm426)
    )
    B::BigInt = zero(BigInt)
    low::BigInt = max(j123, j453, j426, j156)
    high::BigInt = min(jpm123 + j453, jpm132 + j426, jpm231 + j156)
    for x = low:high
        B = -B + binomial(x + 1, j123 + 1) * binomial(jpm123, x - j453) *
                 binomial(jpm132, x - j426) * binomial(jpm231, x - j156)
    end
    return SqrtRational(iphase(high) * B // (dj4 + 1), A)
end

# use 6j-symbol to calculate Racah coefficient
function _dRacah(dj1::BigInt, dj2::BigInt, dj3::BigInt, dj4::BigInt, dj5::BigInt, dj6::BigInt)
    iphase(div(dj1 + dj2 + dj3 + dj4, 2)) * _d6j(dj1, dj2, dj5, dj4, dj3, dj6)
end

# basic 9j-symbol calculation funciton
function _d9j(dj1::BigInt, dj2::BigInt, dj3::BigInt,
    dj4::BigInt, dj5::BigInt, dj6::BigInt,
    dj7::BigInt, dj8::BigInt, dj9::BigInt)
    check_9j(dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9) || return zero(SqrtRational)
    j123::BigInt = div(dj1 + dj2 + dj3, 2)
    j456::BigInt = div(dj4 + dj5 + dj6, 2)
    j789::BigInt = div(dj7 + dj8 + dj9, 2)
    j147::BigInt = div(dj1 + dj4 + dj7, 2)
    j258::BigInt = div(dj2 + dj5 + dj8, 2)
    j369::BigInt = div(dj3 + dj6 + dj9, 2)
    pm123::BigInt = div(dj1 + dj2 - dj3, 2)
    pm132::BigInt = div(dj1 + dj3 - dj2, 2)
    pm231::BigInt = div(dj2 + dj3 - dj1, 2)
    pm456::BigInt = div(dj4 + dj5 - dj6, 2)
    pm465::BigInt = div(dj4 + dj6 - dj5, 2)
    pm564::BigInt = div(dj5 + dj6 - dj4, 2)
    pm789::BigInt = div(dj7 + dj8 - dj9, 2)
    pm798::BigInt = div(dj7 + dj9 - dj8, 2)
    pm897::BigInt = div(dj8 + dj9 - dj7, 2)
    # value in sqrt
    P0_nu::BigInt = binomial(j123 + 1, dj1 + 1) * binomial(dj1, pm123) *
                    binomial(j456 + 1, dj5 + 1) * binomial(dj5, pm456) *
                    binomial(j789 + 1, dj9 + 1) * binomial(dj9, pm798)
    P0_de::BigInt = binomial(j147 + 1, dj1 + 1) * binomial(dj1, div(dj1 + dj4 - dj7, 2)) *
                    binomial(j258 + 1, dj5 + 1) * binomial(dj5, div(dj2 + dj5 - dj8, 2)) *
                    binomial(j369 + 1, dj9 + 1) * binomial(dj9, div(dj3 + dj9 - dj6, 2))
    P0 = P0_nu // P0_de
    dtl::BigInt = max(abs(dj2 - dj6), abs(dj4 - dj8), abs(dj1 - dj9))
    dth::BigInt = min(dj2 + dj6, dj4 + dj8, dj1 + dj9)
    PABC::Rational{BigInt} = zero(Rational{BigInt})
    for dt::BigInt = dtl:2:dth
        j19t::BigInt = div(dj1 + dj9 + dt, 2)
        j26t::BigInt = div(dj2 + dj6 + dt, 2)
        j48t::BigInt = div(dj4 + dj8 + dt, 2)
        Pt_de = binomial(j19t + 1, dt + 1) * binomial(dt, div(dj1 + dt - dj9, 2)) *
                binomial(j26t + 1, dt + 1) * binomial(dt, div(dj2 + dt - dj6, 2)) *
                binomial(j48t + 1, dt + 1) * binomial(dt, div(dj4 + dt - dj8, 2))
        Pt_de *= (dt + 1)^2
        xl::BigInt = max(j123, j369, j26t, j19t)
        xh::BigInt = min(pm123 + j369, pm132 + j26t, pm231 + j19t)
        At::BigInt = zero(BigInt)
        for x = xl:xh
            At = -At + binomial(x + 1, j123 + 1) * binomial(pm123, x - j369) * binomial(pm132, x - j26t) * binomial(pm231, x - j19t)
        end
        yl::BigInt = max(j456, j26t, j258, j48t)
        yh::BigInt = min(pm456 + j26t, pm465 + j258, pm564 + j48t)
        Bt::BigInt = zero(BigInt)
        for y = yl:yh
            Bt = -Bt + binomial(y + 1, j456 + 1) * binomial(pm456, y - j26t) * binomial(pm465, y - j258) * binomial(pm564, y - j48t)
        end
        zl::BigInt = max(j789, j19t, j48t, j147)
        zh::BigInt = min(pm789 + j19t, pm798 + j48t, pm897 + j147)
        Ct::BigInt = zero(BigInt)
        for z = zl:zh
            Ct = -Ct + binomial(z + 1, j789 + 1) * binomial(pm789, z - j19t) * binomial(pm798, z - j48t) * binomial(pm897, z - j147)
        end
        PABC += (iphase(xh + yh + zh) * At * Bt * Ct) // Pt_de
    end
    return SqrtRational(iphase(dth) * PABC, P0)
end

@inline function _lsjj_helper(l1::BigInt, l2::BigInt, dj1::BigInt, dj2::BigInt, J::BigInt)
    pj = div(dj1 + dj2, 2)
    mj = div(dj1 - dj2, 2)
    # j1 = l1 + 1//2, j2 = l2 + 1//2
    (dj1 > 2l1 && dj2 > 2l2) && return exact_sqrt(((pj + J + 1) * (pj - J)) // (2dj1 * dj2))
    # j1 = l1 + 1//2, j2 = l2 - 1//2
    (dj1 > 2l1 && dj2 < 2l2) && return exact_sqrt(((mj + J) * (J - mj + 1)) // (2dj1 * (dj2+2)))
    # j1 = l1 - 1//2, j2 = l2 + 1//2
    (dj1 < 2l1 && dj2 > 2l2) && return -exact_sqrt(((mj + J + 1) * (J - mj)) // (2(dj1+2) * dj2))
    # j1 = l1 - 1//2, j2 = l2 - 1//2
    return exact_sqrt(((pj + J + 2) * (pj - J + 1)) // (2(dj1+2) * (dj2+2)))
end

# S = 0
@inline function _lsjj_S0(l1::BigInt, l2::BigInt, dj1::BigInt, dj2::BigInt, J::BigInt)
    _lsjj_helper(l1, l2, dj1, dj2, J)
end

# S = 1, J = L - 1
function _lsjj_S1_m1(l1::BigInt, l2::BigInt, dj1::BigInt, dj2::BigInt, J::BigInt)
    L = J + 1
    pj = div(dj1 + dj2, 2)
    mj = div(dj1 - dj2, 2)
    pl = l1 + l2
    ml = l1 - l2
    f0 = (J + 1) * (2J + 1)
    fJ = (L + ml) * (L - ml) * (L + pl + 1) * (pl - L + 1)
    fL = (L + mj) * (L - mj) * (L + pj + 1) * (pj - L + 1)
    exact_sqrt(fJ // f0) * _lsjj_helper(l1, l2, dj1, dj2, J) - exact_sqrt(fL // f0) * _lsjj_helper(l1, l2, dj1, dj2, L)
end

# S = 1, J = L
function _lsjj_S1_0(l1::BigInt, l2::BigInt, dj1::BigInt, dj2::BigInt, J::BigInt)
    pj = div(dj1 + dj2, 2)
    mj = div(dj1 - dj2, 2)
    pl = l1 + l2
    ml = l1 - l2
    (mj * (pj + 1) - ml * (pl + 1)) * _lsjj_helper(l1, l2, dj1, dj2, J) / exact_sqrt(J * (J + 1))
end

# S = 1, J = L + 1
function _lsjj_S1_p1(l1::BigInt, l2::BigInt, dj1::BigInt, dj2::BigInt, J::BigInt)
    L = J - 1
    pj = div(dj1 + dj2, 2)
    mj = div(dj1 - dj2, 2)
    pl = l1 + l2
    ml = l1 - l2
    f0 = J * (2J + 1)
    fL = (J + mj) * (J - mj) * (J + pj + 1) * (pj - J + 1)
    fJ = (J + ml) * (J - ml) * (J + pl + 1) * (pl - J + 1)
    exact_sqrt(fL // f0) * _lsjj_helper(l1, l2, dj1, dj2, L) - exact_sqrt(fJ // f0) * _lsjj_helper(l1, l2, dj1, dj2, J)
end

function _lsjj(l1::BigInt, l2::BigInt, dj1::BigInt, dj2::BigInt, L::BigInt, S::BigInt, J::BigInt)
    if abs(dj1 - 2l1) != 1 || abs(dj2 - 2l2) != 1
        return zero(SqrtRational)
    end
    check_couple(2l1, 2l2, 2L) || return zero(SqrtRational)
    check_couple(dj1, dj2, 2J) || return zero(SqrtRational)
    check_couple(2L, 2S, 2J) || return zero(SqrtRational)
    S == 0 && return _lsjj_S0(l1, l2, dj1, dj2, J)
    if S == 1
        J == L - 1 && return _lsjj_S1_m1(l1, l2, dj1, dj2, J)
        J == L && return _lsjj_S1_0(l1, l2, dj1, dj2, J)
        J == L + 1 && return _lsjj_S1_p1(l1, l2, dj1, dj2, J)
    end
    return zero(SqrtRational)
end

function _Moshinsky(N::BigInt, L::BigInt, n::BigInt, l::BigInt, n1::BigInt, l1::BigInt, n2::BigInt, l2::BigInt, Λ::BigInt)
    f1 = 2 * n1 + l1
    f2 = 2 * n2 + l2
    F = 2 * N + L
    f = 2 * n + l
    f1 + f2 == F + f || return zero(SqrtRational)
    χ = f1 + f2
    nl1 = n1 + l1
    nl2 = n2 + l2
    NL = N + L
    nl = n + l
    r1 = binomial(2nl1 + 1, nl1) // (binomial(f1 + 2, n1) * (nl1 + 2))
    r2 = binomial(2nl2 + 1, nl2) // (binomial(f2 + 2, n2) * (nl2 + 2))
    R = binomial(2NL + 1, NL) // (binomial(F + 2, N) * (NL + 2))
    r = binomial(2nl + 1, nl) // (binomial(f + 2, n) * (nl + 2))
    pre_sum = exact_sqrt(r1 * r2 * R * r // big(2)^χ)
    half_lsum = div(l1 + l2 + L + l, 2)
    sum = zero(SqrtRational)
    for fa = 0:min(f1, F)
        fb = f1 - fa
        fc = F - fa
        fd = f2 - fc
        fd >= 0 || continue
        t = exact_sqrt(binomial(f1 + 2, fa + 1) * binomial(f2 + 2, fc + 1) * binomial(F + 2, fa + 1) * binomial(f + 2, fb + 1))
        for la = (fa&0x01):2:fa
            na = div(fa - la, 2)
            nla = na + la
            ta = ((2 * la + 1) * binomial(fa + 1, na)) // binomial(2nla + 1, nla)
            for lb = abs(l1 - la):2:min(la + l1, fb)
                nb = div(fb - lb, 2)
                nlb = nb + lb
                tb = ((2 * lb + 1) * binomial(fb + 1, nb)) // binomial(2nlb + 1, nlb)
                g1 = div(la + lb + l1, 2)
                pCG_ab = SqrtRational(binomial(g1, l1) * binomial(l1, g1 - la), big(1) // (binomial(2g1 + 1, 2(g1 - l1)) * binomial(2l1, 2(g1 - la))))
                for lc = abs(L - la):2:min(la + L, fc)
                    nc = div(fc - lc, 2)
                    nlc = nc + lc
                    tc = ((2 * lc + 1) * binomial(fc + 1, nc)) // binomial(2nlc + 1, nlc)
                    G = div(la + lc + L, 2)
                    pCG_ac = SqrtRational(binomial(G, L) * binomial(L, G - la), big(1) // (binomial(2G + 1, 2(G - L)) * binomial(2L, 2(G - la))))
                    ld_min = max(abs(l2 - lc), abs(l - lb))
                    ld_max = min(fd, lb + l, lc + l2)
                    for ld = ld_min:2:ld_max
                        nd = div(fd - ld, 2)
                        nld = nd + ld
                        td = ((2 * ld + 1) * binomial(fd + 1, nd)) // binomial(2nld + 1, nld)
                        g2 = div(lc + ld + l2, 2)
                        pCG_cd = SqrtRational(binomial(g2, l2) * binomial(l2, g2 - lc), big(1) // (binomial(2g2 + 1, 2(g2 - l2)) * binomial(2l2, 2(g2 - lc))))
                        g = div(lb + ld + l, 2)
                        pCG_bd = SqrtRational(binomial(g, l) * binomial(l, g - lb), big(1) // (binomial(2g + 1, 2(g - l)) * binomial(2l, 2(g - lb))))
                        l_diff = la + lb + lc + ld - half_lsum
                        phase = l_diff >= 0 ? iphase(ld) * 2^l_diff : iphase(ld) // 2^(-l_diff)
                        ninej = d9j(2 * la, 2 * lb, 2 * l1, 2 * lc, 2 * ld, 2 * l2, 2 * L, 2 * l, 2 * Λ)
                        sum = sum + phase * t * ta * tb * tc * td * pCG_ab * pCG_ac * pCG_bd * pCG_cd * ninej
                    end
                end
            end
        end
    end
    return pre_sum * sum
end


"""
    dCG(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer)
CG coefficient function with double angular monentum number parameters, so that the parameters can be integer.
You can calculate `dCG(1, 1, 2, 1, 1, 2)` to calculate the real `CG(1//2, 1//2, 1, 1/2, 1//2, 1)`
"""
@inline dCG(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer) = _dCG(BigInt.((dj1, dj2, dj3, dm1, dm2, dm3))...)

"""
    d3j(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer)
3j-symbol function with double angular monentum number parameters, so that the parameters can be integer.
"""
@inline d3j(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer) = _d3j(BigInt.((dj1, dj2, dj3, dm1, dm2, dm3))...)

"""
    d6j(dj1::Integer, dj2::Integer, dj3::Integer, dj4::Integer, dj5::Integer, dj6::Integer)
6j-symbol function with double angular monentum number parameters, so that the parameters can be integer.
"""
@inline d6j(dj1::Integer, dj2::Integer, dj3::Integer, dj4::Integer, dj5::Integer, dj6::Integer) = _d6j(BigInt.((dj1, dj2, dj3, dj4, dj5, dj6))...)

"""
    dRacah(dj1::Integer, dj2::Integer, dj3::Integer, dj4::Integer, dj5::Integer, dj6::Integer)
Racah coefficient function with double angular momentum parameters, so that the parameters can be integer.
"""
@inline dRacah(dj1::Integer, dj2::Integer, dj3::Integer, dj4::Integer, dj5::Integer, dj6::Integer) = _dRacah(BigInt.((dj1, dj2, dj3, dj4, dj5, dj6))...)

"""
    d9j(dj1::Integer, dj2::Integer, dj3::Integer,
        dj4::Integer, dj5::Integer, dj6::Integer,
        dj7::Integer, dj8::Integer, dj9::Integer)
9j-symbol function with double angular monentum number parameters, so that the parameters can be integer.
"""
@inline d9j(dj1::Integer, dj2::Integer, dj3::Integer,
    dj4::Integer, dj5::Integer, dj6::Integer,
    dj7::Integer, dj8::Integer, dj9::Integer) = _d9j(BigInt.((dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9))...)

@doc raw"""
    CG(j1::Real, j2::Real, j3::Real, m1::Real, m2::Real, m3::Real)
CG coefficient ``\langle j_1m_1 j_2m_2 | j_3m_3 \rangle``
"""
@inline CG(j1::Real, j2::Real, j3::Real, m1::Real, m2::Real, m3::Real) = simplify(_dCG(BigInt.((2j1, 2j2, 2j3, 2m1, 2m2, 2m3))...))

@doc raw"""
    CG0(j1::Integer, j2::Integer, j3::Integer)
CG coefficient special case: ``\langle j_1 0 j_2 0 | j_3 0 \rangle``.
"""
@inline CG0(j1::Integer, j2::Integer, j3::Integer) = simplify(_CG0(big(j1), big(j2), big(j3)))

@doc raw"""
    threeJ(j1::Real, j2::Real, j3::Real, m1::Real, m2::Real, m3::Real)
Wigner 3j-symbol
```math
\begin{pmatrix}
j_1 & j_2 & j_3 \\
m_1 & m_2 & m_3
\end{pmatrix}
```
"""
@inline threeJ(j1::Real, j2::Real, j3::Real, m1::Real, m2::Real, m3::Real) = simplify(_d3j(BigInt.((2j1, 2j2, 2j3, 2m1, 2m2, 2m3))...))

@doc raw"""
    sixJ(j1::Real, j2::Real, j3::Real, j4::Real, j5::Real, j6::Real)
Wigner 6j-symbol
```math
\begin{Bmatrix}
j_1 & j_2 & j_3 \\
j_4 & j_5 & j_6
\end{Bmatrix}
```
"""
@inline sixJ(j1::Real, j2::Real, j3::Real, j4::Real, j5::Real, j6::Real) = simplify(_d6j(BigInt.((2j1, 2j2, 2j3, 2j4, 2j5, 2j6))...))

@doc raw"""
    Racah(j1::Real, j2::Real, j3::Real, j4::Real, j5::Real, j6::Real)
Racah coefficient
```math
W(j_1j_2j_3j_4, j_5j_6) = (-1)^{j_1+j_2+j_3+j_4} \begin{Bmatrix}
j_1 & j_2 & j_5 \\
j_4 & j_3 & j_6
\end{Bmatrix}
```
"""
@inline Racah(j1::Real, j2::Real, j3::Real, j4::Real, j5::Real, j6::Real) = simplify(_dRacah(BigInt.((2j1, 2j2, 2j3, 2j4, 2j5, 2j6))...))

@doc raw"""
    nineJ(j1::Real, j2::Real, j3::Real,
          j4::Real, j5::Real, j6::Real,
          j7::Real, j8::Real, j9::Real)
Wigner 9j-symbol
```math
\begin{Bmatrix}
j_1 & j_2 & j_3 \\
j_4 & j_5 & j_6 \\
j_7 & j_8 & j_9
\end{Bmatrix}
```
"""
@inline nineJ(j1::Real, j2::Real, j3::Real,
    j4::Real, j5::Real, j6::Real,
    j7::Real, j8::Real, j9::Real) = simplify(_d9j(BigInt.((2j1, 2j2, 2j3, 2j4, 2j5, 2j6, 2j7, 2j8, 2j9))...))

@doc raw"""
    norm9J(j1::Real, j2::Real, j3::Real,
           j4::Real, j5::Real, j6::Real,
           j7::Real, j8::Real, j9::Real)
normalized Wigner 9j-symbol
```math
\begin{bmatrix}
j_1 & j_2 & j_3 \\
j_4 & j_5 & j_6 \\
j_7 & j_8 & j_9
\end{bmatrix} = \sqrt{(2j_3+1)(2j_6+1)(2j_7+1)(2j_8+1)} \begin{Bmatrix}
j_1 & j_2 & j_3 \\
j_4 & j_5 & j_6 \\
j_7 & j_8 & j_9
\end{Bmatrix}
```
"""
@inline norm9J(j1::Real, j2::Real, j3::Real,
    j4::Real, j5::Real, j6::Real,
    j7::Real, j8::Real, j9::Real) = simplify(exact_sqrt((2j3 + 1) * (2j6 + 1) * (2j7 + 1) * (2j8 + 1)) * _d9j(BigInt.((2j1, 2j2, 2j3, 2j4, 2j5, 2j6, 2j7, 2j8, 2j9))...))


@doc raw"""
    lsjj(l1::Integer, l2::Integer, j1::Real, j2::Real, L::Integer, S::Integer, J::Integer)
LS-coupling to jj-coupling transformation coefficient
```math
|l_1 l_2 j_1 j_2; J\rangle = \sum_{LS} \langle l_1 l_2 LSJ | l_1 l_2 j_1 j_2; J \rangle |l_1 l_2 LSJ \rangle
```
"""
@inline lsjj(l1::Integer, l2::Integer, j1::Real, j2::Real, L::Integer, S::Integer, J::Integer) = simplify(_lsjj(BigInt.((l1, l2, 2j1, 2j2, L, S, J))...))


@doc raw"""
    Moshinsky(N::Integer, L::Integer, n::Integer, l::Integer, n1::Integer, l1::Integer, n2::Integer, l2::Integer, Λ::Integer)
Moshinsky bracket，Ref: Buck et al. Nuc. Phys. A 600 (1996) 387-402.
Since this function is designed for demonstration the exact result,
It only calculate the case of ``\tan(\beta) = 1``.
"""
@inline Moshinsky(N::Integer, L::Integer, n::Integer, l::Integer, n1::Integer, l1::Integer, n2::Integer, l2::Integer, Λ::Integer) = simplify(_Moshinsky(BigInt.((N, L, n, l, n1, l1, n2, l2, Λ))...))
