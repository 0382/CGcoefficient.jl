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
    for dt::Int = dtl:2:dth
        j19t::Int = div(dj1 + dj9 + dt, 2)
        j26t::Int = div(dj2 + dj6 + dt, 2)
        j48t::Int = div(dj4 + dj8 + dt, 2)
        MPZ.set_ui!(Pt_de, convert(Culong, dt + 1))
        MPZ.mul_ui!(Pt_de, convert(Culong, dt + 1))
        MPZ.mul!(Pt_de, _bigbin(t, j19t + 1, dt + 1))
        MPZ.mul!(Pt_de, _bigbin(t, dt, div(dj1 + dt - dj9, 2)))
        MPZ.mul!(Pt_de, _bigbin(t, j26t + 1, dt + 1))
        MPZ.mul!(Pt_de, _bigbin(t, dt, div(dj2 + dt - dj6, 2)))
        MPZ.mul!(Pt_de, _bigbin(t, j48t + 1, dt + 1))
        MPZ.mul!(Pt_de, _bigbin(t, dt, div(dj4 + dt - dj8, 2)))
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
        MPZ.mul!(At, Bt)
        zl::Int = max(j789, j19t, j48t, j147)
        zh::Int = min(pm789 + j19t, pm798 + j48t, pm897 + j147)
        MPZ.set_ui!(Bt, 0)
        for z = zl:zh
            _bigbin(tx, z + 1, j789 + 1)
            MPZ.mul!(tx, _bigbin(t, pm789, z - j19t))
            MPZ.mul!(tx, _bigbin(t, pm798, z - j48t))
            MPZ.mul!(tx, _bigbin(t, pm897, z - j147))
            MPZ.sub!(Bt, tx, Bt)
        end
        MPZ.mul!(At, Bt)
        if isodd(xh + yh + zh)
            MPZ.neg!(At)
        end
        _divgcd!(tx, At, Pt_de)
        __gmpq_add!(tx, PABC, At, Pt_de)
    end
    if isodd(dth)
        MPZ.neg!(PABC.num)
    end
    Jmax::Int = max(dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9)
    simplify4!(tx, PABC.num, PABC.den, P0_nu, P0_de, cld(5Jmax, 2) + 1)
    return SqrtRational(PABC, Base.unsafe_rational(P0_nu, P0_de))
end

function _omega(ans::BigInt, temp::BigInt, j1::Int, j3::Int, g::Int)
    _bigbin(ans, g, j3)
    MPZ.mul!(ans, _bigbin(temp, j3, g - j1))
end

function _m9j!(
    sum::Rational{BigInt}, temp::BigInt, tx::BigInt, At::BigInt, Bt::BigInt, Pt_de::BigInt,
    j1::Int, j2::Int, j3::Int, j4::Int, j5::Int, j6::Int, j7::Int, j8::Int, j9::Int)
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
    MPZ.set_ui!(sum.num, 0)
    MPZ.set_ui!(sum.den, 1)
    tl::Int = max(abs(j2 - j6), abs(j4 - j8), abs(j1 - j9))
    th::Int = min(j2 + j6, j4 + j8, j1 + j9)
    for t::Int = tl:th
        j19t::Int = j1 + j9 + t
        j26t::Int = j2 + j6 + t
        j48t::Int = j4 + j8 + t
        dt::Int = 2t
        MPZ.set_ui!(Pt_de, convert(Culong, dt + 1))
        MPZ.mul_ui!(Pt_de, convert(Culong, dt + 1))
        MPZ.mul!(Pt_de, _bigbin(temp, j19t + 1, dt + 1))
        MPZ.mul!(Pt_de, _bigbin(temp, dt, j1 + t - j9))
        MPZ.mul!(Pt_de, _bigbin(temp, j26t + 1, dt + 1))
        MPZ.mul!(Pt_de, _bigbin(temp, dt, j2 + t - j6))
        MPZ.mul!(Pt_de, _bigbin(temp, j48t + 1, dt + 1))
        MPZ.mul!(Pt_de, _bigbin(temp, dt, j4 + t - j8))
        xl::Int = max(j123, j369, j26t, j19t)
        xh::Int = min(pm123 + j369, pm132 + j26t, pm231 + j19t)
        MPZ.set_ui!(At, 0)
        for x = xl:xh
            _bigbin(tx, x + 1, j123 + 1)
            MPZ.mul!(tx, _bigbin(temp, pm123, x - j369))
            MPZ.mul!(tx, _bigbin(temp, pm132, x - j26t))
            MPZ.mul!(tx, _bigbin(temp, pm231, x - j19t))
            MPZ.sub!(At, tx, At)
        end
        yl::Int = max(j456, j26t, j258, j48t)
        yh::Int = min(pm456 + j26t, pm465 + j258, pm564 + j48t)
        MPZ.set_ui!(Bt, 0)
        for y = yl:yh
            _bigbin(tx, y + 1, j456 + 1)
            MPZ.mul!(tx, _bigbin(temp, pm456, y - j26t))
            MPZ.mul!(tx, _bigbin(temp, pm465, y - j258))
            MPZ.mul!(tx, _bigbin(temp, pm564, y - j48t))
            MPZ.sub!(Bt, tx, Bt)
        end
        MPZ.mul!(At, Bt)
        zl::Int = max(j789, j19t, j48t, j147)
        zh::Int = min(pm789 + j19t, pm798 + j48t, pm897 + j147)
        MPZ.set_ui!(Bt, 0)
        for z = zl:zh
            _bigbin(tx, z + 1, j789 + 1)
            MPZ.mul!(tx, _bigbin(temp, pm789, z - j19t))
            MPZ.mul!(tx, _bigbin(temp, pm798, z - j48t))
            MPZ.mul!(tx, _bigbin(temp, pm897, z - j147))
            MPZ.sub!(Bt, tx, Bt)
        end
        MPZ.mul!(At, Bt)
        if isodd(xh + yh + zh)
            MPZ.neg!(At)
        end
        _divgcd!(tx, At, Pt_de)
        __gmpq_add!(tx, sum, At, Pt_de)
    end
    return
end

function _lsjj(l1::Int, l2::Int, dj1::Int, dj2::Int, L::Int, S::Int, J::Int)::Rational{Int}
    check_couple(2l1, 2l2, 2L) || return zero(Rational{Int})
    check_couple(dj1, dj2, 2J) || return zero(Rational{Int})
    LSJcase = 0
    if S == 0
        LSJcase = 1
    elseif S == 1
        (L == J - 1) && (LSJcase = 2)
        (J == L) && (LSJcase = 3)
        (L == J + 1) && (LSJcase = 4)
    else
        return zero(Rational{Int})
    end
    m = l1 - l2
    p = l1 + l2
    δ1 = dj1 - 2l1
    δ2 = dj2 - 2l2
    d = 2*(2l1+1)*(2l2+1)
    r = zero(Rational{Int})
    if δ1 == 1 && δ2 == 1
        if LSJcase == 1
            r = (J + p + 2) * (p + 1 - J) // d
        elseif LSJcase == 2
            r = (J + m) * (J - m) // (J * (2J + 1))
            r *= (J + p + 1) * (J + p + 2) // d
        elseif LSJcase == 3
            m == 0 && return zero(Rational{Int})
            r = m * abs(m) // (J * (J + 1))
            r *= (J + p + 2) * (p + 1 - J) // d
        elseif LSJcase == 4
            r = -(L + m) * (L - m) // (L * (2J + 1))
            r *= (p - J) * (p + 1 - J) // d
        end
    elseif δ1 == 1 && δ2 == -1
        if LSJcase == 1
            r = (J + m + 1) * (J - m) // d
        elseif LSJcase == 2
            r = -(p + 1 + J) * (p + 1 - J) // (J * (2J + 1))
            r *= (J + m + 1) * (J + m) // d
        elseif LSJcase == 3
            r = (p + 1) * (p + 1) // (J * (J + 1))
            r *= (J + m + 1) * (J - m) // d
        elseif LSJcase == 4
            r = -(p + 1 + L) * (p + 1 - L) // (L * (2J + 1))
            r *= (J - m + 1) * (J - m) // d
        end
    elseif δ1 == -1 && δ2 == 1
        if LSJcase == 1
            r = -(J + m) * (J - m + 1) // d
        elseif LSJcase == 2
            r = (p + 1 + J) * (p + 1 - J) // (J * (2J + 1))
            r *= (J - m) * (J - m + 1) // d
        elseif LSJcase == 3
            r = (p + 1) * (p + 1) // (J * (J + 1))
            r *= (J + m) * (J - m + 1) // d
        elseif LSJcase == 4
            r = (p + 1 + L) * (p + 1 - L) // (L * (2J + 1))
            r *= (J + m + 1) * (J + m) // d
        end
    elseif δ1 == -1 && δ2 == -1
        if LSJcase == 1
            r = (J + p + 1) * (p - J) // d
        elseif LSJcase == 2
            r = -(J + m) * (J - m) // (J * (2J + 1))
            r *= (p - J) * (p + 1 - J) // d
        elseif LSJcase == 3
            m == 0 && return zero(Rational{Int})
            r = -m * abs(m) // (J * (J + 1))
            r *= (J + p + 1) * (p - J) // d
        elseif LSJcase == 4
            r = (L + m) * (L - m) // (L * (2J + 1))
            r *= (J + p + 1) * (J + p + 2) // d
        end
    end
    return r
end

function _Moshinsky(N::Int, L::Int, n::Int, l::Int, n1::Int, l1::Int, n2::Int, l2::Int, Λ::Int)
    check_couple(2l1, 2l2, 2Λ) || return zero(SqrtRational{BigInt})
    check_couple(2L, 2l, 2Λ) || return zero(SqrtRational{BigInt})
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

    temp = BigInt()
    tx = BigInt()
    At = BigInt()
    Bt = BigInt()
    Pt_de = BigInt()
    nu_a = BigInt()
    de_a = BigInt()
    nu_b = BigInt()
    de_b = BigInt()
    nu_c = BigInt()
    de_c = BigInt()
    nu_d = BigInt()
    de_d = BigInt()
    nu_fa = BigInt()
    M9j = zero(Rational{BigInt})
    
    Snu = BigInt()
    Sde = BigInt()

    _bigbin(Snu, χ + 2, f1 + 1)
    MPZ.mul!(Snu, _bigbin(temp, L + l + Λ + 1, 2Λ + 1))
    MPZ.mul!(Snu, _bigbin(temp, 2Λ, Λ + L - l))
    _bigbin(Sde, χ + 2, F + 1)
    MPZ.mul!(Sde, _bigbin(temp, l1 + l2 + Λ + 1, 2Λ + 1))
    MPZ.mul!(Sde, _bigbin(temp, 2Λ, Λ + l1 - l2))
    _divgcd!(temp, Snu, Sde)

    MPZ.mul_ui!(Snu, convert(Culong, 2l1 + 1))
    MPZ.mul!(Snu, _bigbin(temp, 2nl1 + 1, nl1))
    MPZ.mul!(Sde, _bigbin(temp, f1 + 1, n1))

    MPZ.mul_ui!(Snu, convert(Culong, 2l2 + 1))
    MPZ.mul!(Snu, _bigbin(temp, 2nl2 + 1, nl2))
    MPZ.mul!(Sde, _bigbin(temp, f2 + 1, n2))

    MPZ.mul_ui!(Snu, convert(Culong, 2L + 1))
    MPZ.mul!(Snu, _bigbin(temp, 2NL + 1, NL))
    MPZ.mul!(Sde, _bigbin(temp, F + 1, N))

    MPZ.mul_ui!(Snu, convert(Culong, 2l + 1))
    MPZ.mul!(Snu, _bigbin(temp, 2nl + 1, nl))
    MPZ.mul!(Sde, _bigbin(temp, f + 1, n))

    _divgcd!(temp, Snu, Sde)

    MPZ.set_ui!(tx, convert(Culong, f1 + 2))
    MPZ.mul_ui!(tx, convert(Culong, f2 + 2))
    MPZ.pow_ui!(tx, 2)
    MPZ.mul_2exp!(tx, χ + (l1 + l2 + L + l))    
    MPZ.mul!(Sde, tx)
    _divgcd!(temp, Snu, Sde)

    sum = zero(Rational{BigInt})
    for fa = 0:min(f1, F)
        fb = f1 - fa
        fc = F - fa
        fd = f2 - fc
        fd < 0 && continue
        _bigbin(nu_fa, f1 + 2, fa + 1)
        MPZ.mul!(nu_fa, _bigbin(temp, f2 + 2, fc + 1))
        for la = rem(fa, 2):2:fa
            na = div(fa - la, 2)
            nla = na + la
            MPZ.mul_ui!(nu_a, nu_fa, convert(Culong, 2la + 1))
            MPZ.mul!(nu_a, _bigbin(temp, fa + 1, na))
            MPZ.mul_2exp!(nu_a, la)
            _bigbin(de_a, 2nla + 1, nla)
            _divgcd!(temp, nu_a, de_a)
            for lb = abs(l1 - la):2:min(l1 + la, fb)
                nb = div(fb - lb, 2)
                nlb = nb + lb
                MPZ.mul_ui!(nu_b, nu_a, convert(Culong, 2lb + 1))
                MPZ.mul!(nu_b, _bigbin(temp, fb + 1, nb))
                MPZ.mul_2exp!(nu_b, lb)
                MPZ.mul!(de_b, de_a, _bigbin(temp, 2nlb + 1, nlb))
                _omega(tx, temp, la, l1, div(la + lb + l1, 2))
                MPZ.mul!(nu_b, tx)
                _divgcd!(temp, nu_b, de_b)
                for lc = abs(L - la):2:min(L + la, fc)
                    nc = div(fc - lc, 2)
                    nlc = nc + lc
                    MPZ.mul_ui!(nu_c, nu_b, convert(Culong, 2lc + 1))
                    MPZ.mul!(nu_c, _bigbin(temp, fc + 1, nc))
                    MPZ.mul_2exp!(nu_c, lc)
                    MPZ.mul!(de_c, de_b, _bigbin(temp, 2nlc + 1, nlc))
                    _omega(tx, temp, la, L, div(la + lc + L, 2))
                    MPZ.mul!(nu_c, tx)
                    # Δ(lalcl3)
                    MPZ.set_ui!(tx, convert(Culong, 2L + 1))
                    MPZ.mul!(tx, _bigbin(temp, la + lc + L + 1, 2L + 1))
                    MPZ.mul!(tx, _bigbin(temp, 2L, L + la - lc))
                    MPZ.mul!(de_c, tx)
                    _divgcd!(temp, nu_c, de_c)
                    ldmin = max(abs(l2 - lc), abs(l - lb))
                    ldmax = min(fd, l2 + lc, l + lb)
                    for ld = ldmin:2:ldmax
                        nd = div(fd - ld, 2)
                        nld = nd + ld
                        MPZ.mul_ui!(nu_d, nu_c, convert(Culong, 2ld + 1))
                        MPZ.mul!(nu_d, _bigbin(temp, fd + 1, nd))
                        MPZ.mul_2exp!(nu_d, ld)
                        MPZ.mul!(de_d, de_c, _bigbin(temp, 2nld + 1, nld))
                        _omega(tx, temp, lc, l2, div(lc + ld + l2, 2))
                        MPZ.mul!(nu_d, tx)
                        _omega(tx, temp, lb, l, div(lb + ld + l, 2))
                        MPZ.mul!(nu_d, tx)
                        # Δ(lbldl4)
                        MPZ.set_ui!(tx, convert(Culong, 2l + 1))
                        MPZ.mul!(tx, _bigbin(temp, lb + ld + l + 1, 2l + 1))
                        MPZ.mul!(tx, _bigbin(temp, 2l, l + lb - ld))
                        MPZ.mul!(de_d, tx)
                        _divgcd!(temp, nu_d, de_d)
                        if isodd(ld)
                            MPZ.neg!(nu_d)
                        end
                        _m9j!(M9j, temp, tx, At, Bt, Pt_de, la, lb, l1, lc, ld, l2, L, l, Λ)
                        # don't use MPQ.add! to avoid memory allocation
                        MPZ.mul!(M9j.num, nu_d)
                        MPZ.mul!(M9j.den, de_d)
                        __gmpq_add!(tx, sum, M9j)
                    end
                end
            end
        end
    end
    Jmax = max(f1, f2, F, f, Λ)
    simplify4!(tx, sum.num, sum.den, Snu, Sde, 5Jmax + 1)
    return SqrtRational(sum, Base.unsafe_rational(Snu, Sde))
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
@inline lsjj(l1::Integer, l2::Integer, j1::Real, j2::Real, L::Integer, S::Integer, J::Integer) = begin
    r = _lsjj(Int.((l1, l2, 2j1, 2j2, L, S, J))...)
    return simplify(SqrtRational(sign(r), abs(r)))
end


@doc raw"""
    Moshinsky(N::Integer, L::Integer, n::Integer, l::Integer, n1::Integer, l1::Integer, n2::Integer, l2::Integer, Λ::Integer)
Moshinsky bracket，Ref: Buck et al. Nuc. Phys. A 600 (1996) 387-402.
This function is designed for demonstration the exact result,
It only calculate the case of ``\tan(\beta) = \sqrt{m_1\omega_1/m_2\omega_2} = 1``.
"""
@inline Moshinsky(N::Integer, L::Integer, n::Integer, l::Integer, n1::Integer, l1::Integer, n2::Integer, l2::Integer, Λ::Integer) = _Moshinsky(Int.((N, L, n, l, n1, l1, n2, l2, Λ))...)
