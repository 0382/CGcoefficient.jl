# This file contains the core functions of WignerSymbols and CG-coefficient.

# basic CG coefficient calculation function
function _dCG(dj1::Int, dj2::Int, dj3::Int, dm1::Int, dm2::Int, dm3::Int)
    check_CG(dj1, dj2, dj3, dm1, dm2, dm3) || return zero(SqrtRational{BigInt})
    J::Int = div(dj1 + dj2 + dj3, 2)
    Jm1::Int = J - dj1
    Jm2::Int = J - dj2
    Jm3::Int = J - dj3
    j1mm1::Int = div(dj1 - dm1, 2)
    j2mm2::Int = div(dj2 - dm2, 2)
    j3mm3::Int = div(dj3 - dm3, 2)
    j2pm2::Int = div(dj2 + dm2, 2)
    t = BigInt()
    nA = _bigbin(dj1, Jm2)
    MPZ.mul!(nA, _bigbin(t, dj2, Jm3))
    dA = _bigbin(J + 1, Jm3)
    MPZ.mul!(dA, _bigbin(t, dj1, j1mm1))
    MPZ.mul!(dA, _bigbin(t, dj2, j2mm2))
    MPZ.mul!(dA, _bigbin(t, dj3, j3mm3))

    B::BigInt = zero(BigInt)
    low::Int = max(zero(Int), j1mm1 - Jm2, j2pm2 - Jm1)
    high::Int = min(Jm3, j1mm1, j2pm2)
    tz = BigInt()
    for z in low:high
        _bigbin(tz, Jm3, z)
        MPZ.mul!(tz, _bigbin(t, Jm2, j1mm1 - z))
        MPZ.mul!(tz, _bigbin(t, Jm1, j2pm2 - z))
        MPZ.sub!(B, tz, B)
    end
    if isodd(high)
        MPZ.neg!(B)
    end
    simplify3!(t, B, nA, dA, J + 1)
    return SqrtRational(Base.unsafe_rational(B, t), Base.unsafe_rational(nA, dA))
end

# spaecial case: m1 == m2 == m3 == 0
function _CG0(j1::Int, j2::Int, j3::Int)
    check_couple(2j1, 2j2, 2j3) || return zero(SqrtRational{BigInt})
    J = j1 + j2 + j3
    isodd(J) && return zero(SqrtRational{BigInt})
    g = div(J, 2)
    B = _bigbin(g, j3)
    t = BigInt()
    MPZ.mul!(B, _bigbin(t, j3, g - j1))
    if isodd(g - j3)
        MPZ.neg!(B)
    end
    dA = _bigbin(J + 1, 2j3 + 1)
    MPZ.mul!(dA, _bigbin(t, 2j3, J - 2j1))
    nA = _divgcd!(BigInt(), B, dA)
    _divgcd!(t, nA, dA)
    MPZ.set_ui!(t, 1)
    simplify3!(t, B, nA, dA, J + 1)
    return SqrtRational(Base.unsafe_rational(B, t), Base.unsafe_rational(nA, dA))
end

# basic 3j-symbol calculation funciton
function _d3j(dj1::Int, dj2::Int, dj3::Int, dm1::Int, dm2::Int, dm3::Int)
    check_3j(dj1, dj2, dj3, dm1, dm2, dm3) || return zero(SqrtRational{BigInt})
    J::Int = div(dj1 + dj2 + dj3, 2)
    Jm1::Int = J - dj1
    Jm2::Int = J - dj2
    Jm3::Int = J - dj3
    j1mm1::Int = div(dj1 - dm1, 2)
    j2mm2::Int = div(dj2 - dm2, 2)
    j3mm3::Int = div(dj3 - dm3, 2)
    j1pm1::Int = div(dj1 + dm1, 2)
    t = BigInt()
    nA = _bigbin(dj1, Jm2)
    MPZ.mul!(nA, _bigbin(t, dj2, Jm1))
    dA = _bigbin(J, Jm3)
    MPZ.mul!(dA, _bigbin(t, dj1, j1mm1))
    MPZ.mul!(dA, _bigbin(t, dj2, j2mm2))
    MPZ.mul!(dA, _bigbin(t, dj3, j3mm3))
    MPZ.mul_ui!(dA, convert(Culong, J + 1))

    B::BigInt = zero(BigInt)
    low::Int = max(zero(Int), j1pm1 - Jm2, j2mm2 - Jm1)
    high::Int = min(Jm3, j1pm1, j2mm2)
    tz = BigInt()
    for z in low:high
        _bigbin(tz, Jm3, z)
        MPZ.mul!(tz, _bigbin(t, Jm2, j1pm1 - z))
        MPZ.mul!(tz, _bigbin(t, Jm1, j2mm2 - z))
        MPZ.sub!(B, tz, B)
    end
    if isodd(dj1 + div(dj3 + dm3, 2) + high)
        MPZ.neg!(B)
    end
    simplify3!(t, B, nA, dA, J + 1)
    return SqrtRational(Base.unsafe_rational(B, t), Base.unsafe_rational(nA, dA))
end

# basic 6j-symbol calculation funciton
function _d6j(dj1::Int, dj2::Int, dj3::Int, dj4::Int, dj5::Int, dj6::Int)
    check_6j(dj1, dj2, dj3, dj4, dj5, dj6) || return zero(SqrtRational{BigInt})
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
    # value in sqrt
    t = BigInt()
    nA = _bigbin(j123 + 1, dj1 + 1)
    MPZ.mul!(nA, _bigbin(t, dj1, jpm123))
    dA = _bigbin(j156 + 1, dj1 + 1)
    MPZ.mul!(dA, _bigbin(t, dj1, jpm156))
    MPZ.mul!(dA, _bigbin(t, j453 + 1, dj4 + 1))
    MPZ.mul!(dA, _bigbin(t, dj4, jpm453))
    MPZ.mul!(dA, _bigbin(t, j426 + 1, dj4 + 1))
    MPZ.mul!(dA, _bigbin(t, dj4, jpm426))
    MPZ.mul_ui!(dA, convert(Culong, dj4 + 1))
    MPZ.mul_ui!(dA, convert(Culong, dj4 + 1))
    B::BigInt = zero(BigInt)
    low::Int = max(j123, j453, j426, j156)
    high::Int = min(jpm123 + j453, jpm132 + j426, jpm231 + j156)
    tx = BigInt()
    for x = low:high
        _bigbin(t, x + 1, j123 + 1)
        MPZ.mul!(t, _bigbin(tx, jpm123, x - j453))
        MPZ.mul!(t, _bigbin(tx, jpm132, x - j426))
        MPZ.mul!(t, _bigbin(tx, jpm231, x - j156))
        MPZ.sub!(B, t, B)
    end
    if isodd(high)
        MPZ.neg!(B)
    end
    simplify3!(t, B, nA, dA, high + 1)
    return SqrtRational(Base.unsafe_rational(B, t), Base.unsafe_rational(nA, dA))
end

# use 6j-symbol to calculate Racah coefficient
function _dRacah(dj1::Int, dj2::Int, dj3::Int, dj4::Int, dj5::Int, dj6::Int)
    ans = _d6j(dj1, dj2, dj5, dj4, dj3, dj6)
    if isodd(div(dj1 + dj2 + dj3 + dj4, 2))
        MPZ.neg!(ans.s.num)
    end
    return ans
end

# basic 9j-symbol calculation funciton
function _d9j(dj1::Int, dj2::Int, dj3::Int,
    dj4::Int, dj5::Int, dj6::Int,
    dj7::Int, dj8::Int, dj9::Int)
    check_9j(dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9) || return zero(SqrtRational{BigInt})
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
    t = BigInt()
    P0_nu::BigInt = _bigbin(j123 + 1, dj1 + 1)
    MPZ.mul!(P0_nu, _bigbin(t, dj1, pm123))
    MPZ.mul!(P0_nu, _bigbin(t, j456 + 1, dj5 + 1))
    MPZ.mul!(P0_nu, _bigbin(t, dj5, pm456))
    MPZ.mul!(P0_nu, _bigbin(t, j789 + 1, dj9 + 1))
    MPZ.mul!(P0_nu, _bigbin(t, dj9, pm798))
    P0_de::BigInt = _bigbin(j147 + 1, dj1 + 1)
    MPZ.mul!(P0_de, _bigbin(t, dj1, div(dj1 + dj4 - dj7, 2)))
    MPZ.mul!(P0_de, _bigbin(t, j258 + 1, dj5 + 1))
    MPZ.mul!(P0_de, _bigbin(t, dj5, div(dj2 + dj5 - dj8, 2)))
    MPZ.mul!(P0_de, _bigbin(t, j369 + 1, dj9 + 1))
    MPZ.mul!(P0_de, _bigbin(t, dj9, div(dj3 + dj9 - dj6, 2)))
    dtl::Int = max(abs(dj2 - dj6), abs(dj4 - dj8), abs(dj1 - dj9))
    dth::Int = min(dj2 + dj6, dj4 + dj8, dj1 + dj9)
    PABC::Rational{BigInt} = zero(Rational{BigInt})
    Pt_de::BigInt = BigInt()
    tx = BigInt()
    At = BigInt()
    Bt = BigInt()
    Ct = BigInt()
    for dt::Int = dtl:2:dth
        j19t::Int = div(dj1 + dj9 + dt, 2)
        j26t::Int = div(dj2 + dj6 + dt, 2)
        j48t::Int = div(dj4 + dj8 + dt, 2)
        _bigbin(Pt_de, j19t + 1, dt + 1)
        MPZ.mul!(Pt_de, _bigbin(t, dt, div(dj1 + dt - dj9, 2)))
        MPZ.mul!(Pt_de, _bigbin(t, j26t + 1, dt + 1))
        MPZ.mul!(Pt_de, _bigbin(t, dt, div(dj2 + dt - dj6, 2)))
        MPZ.mul!(Pt_de, _bigbin(t, j48t + 1, dt + 1))
        MPZ.mul!(Pt_de, _bigbin(t, dt, div(dj4 + dt - dj8, 2)))
        MPZ.mul_ui!(Pt_de, convert(Culong, dt + 1))
        MPZ.mul_ui!(Pt_de, convert(Culong, dt + 1))
        xl::Int = max(j123, j369, j26t, j19t)
        xh::Int = min(pm123 + j369, pm132 + j26t, pm231 + j19t)
        MPZ.set_ui!(At, 0)
        for x = xl:xh
            _bigbin(tx, x + 1, j123 + 1)
            MPZ.mul!(tx, _bigbin(t, pm123, x - j369))
            MPZ.mul!(tx, _bigbin(t, pm132, x - j26t))
            MPZ.mul!(tx, _bigbin(t, pm231, x - j19t))
            MPZ.sub!(At, tx, At)
        end
        yl::Int = max(j456, j26t, j258, j48t)
        yh::Int = min(pm456 + j26t, pm465 + j258, pm564 + j48t)
        MPZ.set_ui!(Bt, 0)
        for y = yl:yh
            _bigbin(tx, y + 1, j456 + 1)
            MPZ.mul!(tx, _bigbin(t, pm456, y - j26t))
            MPZ.mul!(tx, _bigbin(t, pm465, y - j258))
            MPZ.mul!(tx, _bigbin(t, pm564, y - j48t))
            MPZ.sub!(Bt, tx, Bt)
        end
        zl::Int = max(j789, j19t, j48t, j147)
        zh::Int = min(pm789 + j19t, pm798 + j48t, pm897 + j147)
        MPZ.set_ui!(Ct, 0)
        for z = zl:zh
            _bigbin(tx, z + 1, j789 + 1)
            MPZ.mul!(tx, _bigbin(t, pm789, z - j19t))
            MPZ.mul!(tx, _bigbin(t, pm798, z - j48t))
            MPZ.mul!(tx, _bigbin(t, pm897, z - j147))
            MPZ.sub!(Ct, tx, Ct)
        end
        MPZ.set!(t, At)
        MPZ.mul!(t, Bt)
        MPZ.mul!(t, Ct)
        if isodd(xh + yh + zh)
            MPZ.neg!(t)
        end
        _divgcd!(tx, t, Pt_de)
        Base.GMP.MPQ.add!(PABC, Base.unsafe_rational(t, Pt_de))
    end
    if isodd(dth)
        MPZ.neg!(PABC.num)
    end
    Jmax::Int = max(dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9)
    simplify4!(t, PABC.num, PABC.den, P0_nu, P0_de, cld(5Jmax, 2) + 1)
    return SqrtRational(PABC, Base.unsafe_rational(P0_nu, P0_de))
end

@inline function _lsjj_helper(l1::Int, l2::Int, dj1::Int, dj2::Int, J::Int)
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
@inline function _lsjj_S0(l1::Int, l2::Int, dj1::Int, dj2::Int, J::Int)
    _lsjj_helper(l1, l2, dj1, dj2, J)
end

# S = 1, J = L - 1
function _lsjj_S1_m1(l1::Int, l2::Int, dj1::Int, dj2::Int, J::Int)
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
function _lsjj_S1_0(l1::Int, l2::Int, dj1::Int, dj2::Int, J::Int)
    pj = div(dj1 + dj2, 2)
    mj = div(dj1 - dj2, 2)
    pl = l1 + l2
    ml = l1 - l2
    _lsjj_helper(l1, l2, dj1, dj2, J) * SqrtRational(mj * (pj + 1) - ml * (pl + 1), 1//(J * (J + 1)))
end

# S = 1, J = L + 1
function _lsjj_S1_p1(l1::Int, l2::Int, dj1::Int, dj2::Int, J::Int)
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

function _lsjj(l1::Int, l2::Int, dj1::Int, dj2::Int, L::Int, S::Int, J::Int)
    if abs(dj1 - 2l1) != 1 || abs(dj2 - 2l2) != 1
        return zero(SqrtRational{Int})
    end
    check_couple(2l1, 2l2, 2L) || return zero(SqrtRational{Int})
    check_couple(dj1, dj2, 2J) || return zero(SqrtRational{Int})
    check_couple(2L, 2S, 2J) || return zero(SqrtRational{Int})
    S == 0 && return _lsjj_S0(l1, l2, dj1, dj2, J)
    if S == 1
        J == L - 1 && return _lsjj_S1_m1(l1, l2, dj1, dj2, J)
        J == L && return _lsjj_S1_0(l1, l2, dj1, dj2, J)
        J == L + 1 && return _lsjj_S1_p1(l1, l2, dj1, dj2, J)
    end
    return zero(SqrtRational{Int})
end

function _Moshinsky(N::Int, L::Int, n::Int, l::Int, n1::Int, l1::Int, n2::Int, l2::Int, Λ::Int)
    f1 = 2 * n1 + l1
    f2 = 2 * n2 + l2
    F = 2 * N + L
    f = 2 * n + l
    f1 + f2 == F + f || return zero(SqrtRational{BigInt})
    χ = f1 + f2
    nl1 = n1 + l1
    nl2 = n2 + l2
    NL = N + L
    nl = n + l
    r1 = _bigbin(2nl1 + 1, nl1) // (_bigbin(f1 + 2, n1) * (nl1 + 2))
    r2 = _bigbin(2nl2 + 1, nl2) // (_bigbin(f2 + 2, n2) * (nl2 + 2))
    R = _bigbin(2NL + 1, NL) // (_bigbin(F + 2, N) * (NL + 2))
    r = _bigbin(2nl + 1, nl) // (_bigbin(f + 2, n) * (nl + 2))
    pre_sum = exact_sqrt(r1 * r2 * R * r // big(2)^χ)
    half_lsum = div(l1 + l2 + L + l, 2)
    sum = zero(SqrtRational{BigInt})
    for fa = 0:min(f1, F)
        fb = f1 - fa
        fc = F - fa
        fd = f2 - fc
        fd >= 0 || continue
        t = exact_sqrt(_bigbin(f1 + 2, fa + 1) * _bigbin(f2 + 2, fc + 1) * _bigbin(F + 2, fa + 1) * _bigbin(f + 2, fb + 1))
        for la = (fa&0x01):2:fa
            na = div(fa - la, 2)
            nla = na + la
            ta = ((2 * la + 1) * _bigbin(fa + 1, na)) // _bigbin(2nla + 1, nla)
            for lb = abs(l1 - la):2:min(la + l1, fb)
                nb = div(fb - lb, 2)
                nlb = nb + lb
                tb = ((2 * lb + 1) * _bigbin(fb + 1, nb)) // _bigbin(2nlb + 1, nlb)
                g1 = div(la + lb + l1, 2)
                pCG_ab = SqrtRational(_bigbin(g1, l1) * _bigbin(l1, g1 - la), big(1) // (_bigbin(2g1 + 1, 2(g1 - l1)) * _bigbin(2l1, 2(g1 - la))))
                for lc = abs(L - la):2:min(la + L, fc)
                    nc = div(fc - lc, 2)
                    nlc = nc + lc
                    tc = ((2 * lc + 1) * _bigbin(fc + 1, nc)) // _bigbin(2nlc + 1, nlc)
                    G = div(la + lc + L, 2)
                    pCG_ac = SqrtRational(_bigbin(G, L) * _bigbin(L, G - la), big(1) // (_bigbin(2G + 1, 2(G - L)) * _bigbin(2L, 2(G - la))))
                    ld_min = max(abs(l2 - lc), abs(l - lb))
                    ld_max = min(fd, lb + l, lc + l2)
                    for ld = ld_min:2:ld_max
                        nd = div(fd - ld, 2)
                        nld = nd + ld
                        td = ((2 * ld + 1) * _bigbin(fd + 1, nd)) // _bigbin(2nld + 1, nld)
                        g2 = div(lc + ld + l2, 2)
                        pCG_cd = SqrtRational(_bigbin(g2, l2) * _bigbin(l2, g2 - lc), big(1) // (_bigbin(2g2 + 1, 2(g2 - l2)) * _bigbin(2l2, 2(g2 - lc))))
                        g = div(lb + ld + l, 2)
                        pCG_bd = SqrtRational(_bigbin(g, l) * _bigbin(l, g - lb), big(1) // (_bigbin(2g + 1, 2(g - l)) * _bigbin(2l, 2(g - lb))))
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
@inline dCG(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer) = _dCG(Int.((dj1, dj2, dj3, dm1, dm2, dm3))...)

"""
    d3j(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer)
3j-symbol function with double angular monentum number parameters, so that the parameters can be integer.
"""
@inline d3j(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer) = _d3j(Int.((dj1, dj2, dj3, dm1, dm2, dm3))...)

"""
    d6j(dj1::Integer, dj2::Integer, dj3::Integer, dj4::Integer, dj5::Integer, dj6::Integer)
6j-symbol function with double angular monentum number parameters, so that the parameters can be integer.
"""
@inline d6j(dj1::Integer, dj2::Integer, dj3::Integer, dj4::Integer, dj5::Integer, dj6::Integer) = _d6j(Int.((dj1, dj2, dj3, dj4, dj5, dj6))...)

"""
    dRacah(dj1::Integer, dj2::Integer, dj3::Integer, dj4::Integer, dj5::Integer, dj6::Integer)
Racah coefficient function with double angular momentum parameters, so that the parameters can be integer.
"""
@inline dRacah(dj1::Integer, dj2::Integer, dj3::Integer, dj4::Integer, dj5::Integer, dj6::Integer) = _dRacah(Int.((dj1, dj2, dj3, dj4, dj5, dj6))...)

"""
    d9j(dj1::Integer, dj2::Integer, dj3::Integer,
        dj4::Integer, dj5::Integer, dj6::Integer,
        dj7::Integer, dj8::Integer, dj9::Integer)
9j-symbol function with double angular monentum number parameters, so that the parameters can be integer.
"""
@inline d9j(dj1::Integer, dj2::Integer, dj3::Integer,
    dj4::Integer, dj5::Integer, dj6::Integer,
    dj7::Integer, dj8::Integer, dj9::Integer) = _d9j(Int.((dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9))...)

@doc raw"""
    CG(j1::Real, j2::Real, j3::Real, m1::Real, m2::Real, m3::Real)
CG coefficient ``\langle j_1m_1 j_2m_2 | j_3m_3 \rangle``
"""
@inline CG(j1::Real, j2::Real, j3::Real, m1::Real, m2::Real, m3::Real) = _dCG(Int.((2j1, 2j2, 2j3, 2m1, 2m2, 2m3))...)

@doc raw"""
    CG0(j1::Integer, j2::Integer, j3::Integer)
CG coefficient special case: ``\langle j_1 0 j_2 0 | j_3 0 \rangle``.
"""
@inline CG0(j1::Integer, j2::Integer, j3::Integer) = _CG0(Int.((j1, j2, j3))...)

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
@inline threeJ(j1::Real, j2::Real, j3::Real, m1::Real, m2::Real, m3::Real) = _d3j(Int.((2j1, 2j2, 2j3, 2m1, 2m2, 2m3))...)

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
@inline sixJ(j1::Real, j2::Real, j3::Real, j4::Real, j5::Real, j6::Real) = _d6j(Int.((2j1, 2j2, 2j3, 2j4, 2j5, 2j6))...)

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
@inline Racah(j1::Real, j2::Real, j3::Real, j4::Real, j5::Real, j6::Real) = _dRacah(Int.((2j1, 2j2, 2j3, 2j4, 2j5, 2j6))...)

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
@inline nineJ(j1::Real, j2::Real, j3::Real, j4::Real, j5::Real, j6::Real, j7::Real, j8::Real, j9::Real) = _d9j(Int.((2j1, 2j2, 2j3, 2j4, 2j5, 2j6, 2j7, 2j8, 2j9))...)

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
@inline norm9J(j1::Real, j2::Real, j3::Real, j4::Real, j5::Real, j6::Real, j7::Real, j8::Real, j9::Real) = begin
    dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9 = Int.((2j1, 2j2, 2j3, 2j4, 2j5, 2j6, 2j7, 2j8, 2j9))
    Jmax = max(dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9)
    ans = _d9j(dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9)
    MPZ.mul_ui!(ans.r.num, convert(Culong, dj3 + 1))
    MPZ.mul_ui!(ans.r.num, convert(Culong, dj6 + 1))
    MPZ.mul_ui!(ans.r.num, convert(Culong, dj7 + 1))
    MPZ.mul_ui!(ans.r.num, convert(Culong, dj8 + 1))
    simplify4!(BigInt(), ans.s.num, ans.s.den, ans.r.num, ans.r.den, cld(5Jmax, 2) + 1)
    return ans
end

@doc raw"""
    lsjj(l1::Integer, l2::Integer, j1::Real, j2::Real, L::Integer, S::Integer, J::Integer)
LS-coupling to jj-coupling transformation coefficient
```math
\langle l_1l_2L,s_1s_2S;J|l_1s_1j_1,l_2s_2j_2;J\rangle = \begin{bmatrix}l_1 & s_1 & j_1 \\ l_2 & s_2 & j_2 \\ L & S & J\end{bmatrix}
```
"""
@inline lsjj(l1::Integer, l2::Integer, j1::Real, j2::Real, L::Integer, S::Integer, J::Integer) = simplify(_lsjj(Int.((l1, l2, 2j1, 2j2, L, S, J))...))


@doc raw"""
    Moshinsky(N::Integer, L::Integer, n::Integer, l::Integer, n1::Integer, l1::Integer, n2::Integer, l2::Integer, Λ::Integer)
Moshinsky bracket，Ref: Buck et al. Nuc. Phys. A 600 (1996) 387-402.
Since this function is designed for demonstration the exact result,
It only calculate the case of ``\tan(\beta) = 1``.
"""
@inline Moshinsky(N::Integer, L::Integer, n::Integer, l::Integer, n1::Integer, l1::Integer, n2::Integer, l2::Integer, Λ::Integer) = simplify(_Moshinsky(Int.((N, L, n, l, n1, l1, n2, l2, Λ))...))
