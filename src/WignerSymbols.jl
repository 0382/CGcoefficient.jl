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
    check_couple(j1, j2, j3) || return zero(SqrtRational{BigInt})
    J = j1 + j2 + j3
    g = div(J, 2)
    sn = _bigbin(g, j3)
    sd = _bigbin(j3, g - j1)
    MPZ.mul!(sn, sd)
    if isodd(g - j3)
        MPZ.neg!(sn)
    end
    rd = _bigbin(J + 1, 2j3 + 1)
    rn = _bigbin(2j3, J - 2j1)
    MPZ.mul!(rd, rn)
    simplify2!(sn, sd, rn, rd, J + 1)
    return SqrtRational(Base.unsafe_rational(sn, sd), Base.unsafe_rational(rn, rd))
end

function _3j0(j1::Int, j2::Int, j3::Int)
    check_couple_(j1, j2, j3) || return zero(SqrtRational{BigInt})
    J = j1 + j2 + j3
    g = div(J, 2)
    sn = _bigbin(g, j3)
    sd = _bigbin(j3, g - j1)
    MPZ.mul!(sn, sd)
    if isodd(g)
        MPZ.neg!(sn)
    end
    rd = _bigbin(J + 1, 2j3 + 1)
    rn = _bigbin(2j3, J - 2j1)
    MPZ.mul!(rd, rn)
    MPZ.mul_ui!(rd, convert(Culong, 2j3 + 1))
    simplify2!(sn, sd, rn, rd, J + 1)
    return SqrtRational(Base.unsafe_rational(sn, sd), Base.unsafe_rational(rn, rd))
end

function _xGaunt(l1::Int, l2::Int, l3::Int, m1::Int, m2::Int, m3::Int)
    check_Gaunt(l1, l2, l3, m1, m2, m3) || return zero(SqrtRational{BigInt})
    J::Int = l1 + l2 + l3
    g::Int = div(J, 2)
    Jm1::Int = l2 + l3 - l1
    Jm2::Int = l1 + l3 - l2
    Jm3::Int = l1 + l2 - l3
    j1mm1::Int = l1 - m1
    j2mm2::Int = l2 - m2
    j3mm3::Int = l3 - m3
    j1pm1::Int = l1 + m1

    sn = zero(BigInt)
    sd = zero(BigInt)
    temp = zero(BigInt)
    low = max(zero(Int), j1pm1 - Jm2, j2mm2 - Jm1)
    high = min(Jm3, j1pm1, j2mm2)
    for z in low:high
        _bigbin(sd, Jm3, z)
        MPZ.mul!(sd, _bigbin(temp, Jm2, j1pm1 - z))
        MPZ.mul!(sd, _bigbin(temp, Jm1, j2mm2 - z))
        MPZ.sub!(sn, sd, sn)
    end
    if iszero(sn)
        return zero(SqrtRational{BigInt})
    end
    if isodd(g + j3mm3 + high)
        MPZ.neg!(sn)
    end
    MPZ.mul!(sn, _bigbin(temp, g, g - l3))
    MPZ.mul!(sn, _bigbin(temp, l3, g - l1))
    rn = zero(BigInt)
    rd = zero(BigInt)
    MPZ.set_ui!(rd, 4)
    MPZ.mul!(rd, _bigbin(temp, 2l1, j1mm1))
    MPZ.mul!(rd, _bigbin(temp, 2l2, j2mm2))
    MPZ.mul!(rd, _bigbin(temp, 2l3, j3mm3))
    MPZ.set_ui!(rn, convert(Culong, 2l3 + 1))
    MPZ.mul!(rn, _bigbin(temp, 2l3, Jm1))
    MPZ.mul!(rd, _bigbin(temp, J + 1, 2l1 + 1))
    MPZ.mul!(rd, _bigbin(temp, J + 1, 2l2 + 1))
    simplify3!(sd, sn, rn, rd, J + 1)
    return SqrtRational(Base.unsafe_rational(sn, sd), Base.unsafe_rational(rn, rd))
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
function _d9j(dj1::Int, dj2::Int, dj3::Int, dj4::Int, dj5::Int, dj6::Int, dj7::Int, dj8::Int, dj9::Int)
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
    dtl = max(abs(dj2 - dj6), abs(dj4 - dj8), abs(dj1 - dj9))
    dth = min(dj2 + dj6, dj4 + dj8, dj1 + dj9)
    sum = zero(BigInt)
    ABC = zero(BigInt)
    Pt = zero(BigInt)
    tx = zero(BigInt)
    temp = zero(BigInt)
    for dt = dtl:2:dth
        j19t = div(dj1 + dj9 + dt, 2)
        j26t = div(dj2 + dj6 + dt, 2)
        j48t = div(dj4 + dj8 + dt, 2)
        MPZ.set_ui!(ABC, convert(Culong, dt + 1))
        xl = max(j123, j369, j26t, j19t)
        xh = min(pm123 + j369, pm132 + j26t, pm231 + j19t)
        MPZ.set_ui!(Pt, 0)
        for x = xl:xh
            _bigbin(tx, x + 1, j19t + 1)
            MPZ.mul!(tx, _bigbin(temp, div(dj1 + dj9 - dt, 2), x - j26t))
            MPZ.mul!(tx, _bigbin(temp, div(dj1 + dt - dj9, 2), x - j369))
            MPZ.mul!(tx, _bigbin(temp, div(dt + dj9 - dj1, 2), x - j123))
            MPZ.sub!(Pt, tx, Pt)
        end
        MPZ.mul!(ABC, Pt)
        yl = max(j456, j26t, j258, j48t)
        yh = min(pm456 + j26t, pm465 + j258, pm564 + j48t)
        MPZ.set_ui!(Pt, 0)
        for y = yl:yh
            _bigbin(tx, y + 1, j26t + 1)
            MPZ.mul!(tx, _bigbin(temp, div(dj2 + dj6 - dt, 2), y - j48t))
            MPZ.mul!(tx, _bigbin(temp, div(dj6 + dt - dj2, 2), y - j258))
            MPZ.mul!(tx, _bigbin(temp, div(dt + dj2 - dj6, 2), y - j456))
            MPZ.sub!(Pt, tx, Pt)
        end
        MPZ.mul!(ABC, Pt)
        zl = max(j789, j19t, j48t, j147)
        zh = min(pm789 + j19t, pm798 + j48t, pm897 + j147)
        MPZ.set_ui!(Pt, 0)
        for z = zl:zh
            _bigbin(tx, z + 1, j48t + 1)
            MPZ.mul!(tx, _bigbin(temp, div(dj4 + dj8 - dt, 2), z - j19t))
            MPZ.mul!(tx, _bigbin(temp, div(dj8 + dt - dj4, 2), z - j147))
            MPZ.mul!(tx, _bigbin(temp, div(dt + dj4 - dj8, 2), z - j789))
            MPZ.sub!(Pt, tx, Pt)
        end
        MPZ.mul!(ABC, Pt)
        if isodd(xh + yh + zh)
            MPZ.neg!(ABC)
        end
        MPZ.add!(sum, ABC)
    end
    if iszero(sum)
        return zero(SqrtRational{BigInt})
    end
    if isodd(dth)
        MPZ.neg!(sum)
    end
    MPZ.set_ui!(Pt, convert(Culong, dj3 + 1))
    MPZ.mul_ui!(Pt, convert(Culong, dj6 + 1))
    MPZ.mul_ui!(Pt, convert(Culong, dj9 + 1))
    MPZ.mul_ui!(Pt, convert(Culong, dj7 + 1))
    MPZ.mul_ui!(Pt, convert(Culong, dj8 + 1))
    MPZ.mul_ui!(Pt, convert(Culong, dj9 + 1))
    MPZ.mul!(Pt, _bigbin(temp, j123 + 1, dj3 + 1))
    MPZ.mul!(Pt, _bigbin(temp, j456 + 1, dj6 + 1))
    MPZ.mul!(Pt, _bigbin(temp, j789 + 1, dj9 + 1))
    MPZ.mul!(Pt, _bigbin(temp, dj3, pm231))
    MPZ.mul!(Pt, _bigbin(temp, dj6, pm564))
    MPZ.mul!(Pt, _bigbin(temp, dj9, pm897))
    MPZ.mul!(Pt, _bigbin(temp, j147 + 1, dj7 + 1))
    MPZ.mul!(Pt, _bigbin(temp, j258 + 1, dj8 + 1))
    MPZ.mul!(Pt, _bigbin(temp, j369 + 1, dj9 + 1))
    MPZ.mul!(Pt, _bigbin(temp, dj7, div(dj7 - dj1 + dj4, 2)))
    MPZ.mul!(Pt, _bigbin(temp, dj8, div(dj8 - dj2 + dj5, 2)))
    MPZ.mul!(Pt, _bigbin(temp, dj9, div(dj9 - dj3 + dj6, 2)))
    Jmax = max(j123, j456, j789, j147, j258, j369)
    simplify2!(sum, tx, ABC, Pt, Jmax + 1)
    return SqrtRational(Base.unsafe_rational(sum, tx), Base.unsafe_rational(ABC, Pt))
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
    end
    iszero(LSJcase) && return zero(Rational{Int})
    m = l1 - l2
    p = l1 + l2
    δ1 = dj1 - 2l1
    δ2 = dj2 - 2l2
    d = 2 * (2l1 + 1) * (2l2 + 1)
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


function _omega(ans::BigInt, temp::BigInt, j1::Int, j3::Int, g::Int)
    _bigbin(ans, g, j3)
    return MPZ.mul!(ans, _bigbin(temp, j3, g - j1))
end

function _m9j!(
    sum::BigInt, temp::BigInt, tx::BigInt, Pt::BigInt, ABC::BigInt,
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
    MPZ.set_ui!(sum, 0)
    tl::Int = max(abs(j2 - j6), abs(j4 - j8), abs(j1 - j9))
    th::Int = min(j2 + j6, j4 + j8, j1 + j9)
    for t::Int = tl:th
        j19t::Int = j1 + j9 + t
        j26t::Int = j2 + j6 + t
        j48t::Int = j4 + j8 + t
        MPZ.set_ui!(ABC, convert(Culong, 2t + 1))
        xl::Int = max(j123, j369, j26t, j19t)
        xh::Int = min(pm123 + j369, pm132 + j26t, pm231 + j19t)
        MPZ.set_ui!(Pt, 0)
        for x = xl:xh
            _bigbin(tx, x + 1, j19t + 1)
            MPZ.mul!(tx, _bigbin(temp, j1 + j9 - t, x - j26t))
            MPZ.mul!(tx, _bigbin(temp, j1 + t - j9, x - j369))
            MPZ.mul!(tx, _bigbin(temp, t + j9 - j1, x - j123))
            MPZ.sub!(Pt, tx, Pt)
        end
        MPZ.mul!(ABC, Pt)
        yl::Int = max(j456, j26t, j258, j48t)
        yh::Int = min(pm456 + j26t, pm465 + j258, pm564 + j48t)
        MPZ.set_ui!(Pt, 0)
        for y = yl:yh
            _bigbin(tx, y + 1, j26t + 1)
            MPZ.mul!(tx, _bigbin(temp, j2 + j6 - t, y - j48t))
            MPZ.mul!(tx, _bigbin(temp, j6 + t - j2, y - j258))
            MPZ.mul!(tx, _bigbin(temp, t + j2 - j6, y - j456))
            MPZ.sub!(Pt, tx, Pt)
        end
        MPZ.mul!(ABC, Pt)
        zl::Int = max(j789, j19t, j48t, j147)
        zh::Int = min(pm789 + j19t, pm798 + j48t, pm897 + j147)
        MPZ.set_ui!(Pt, 0)
        for z = zl:zh
            _bigbin(tx, z + 1, j48t + 1)
            MPZ.mul!(tx, _bigbin(temp, j4 + j8 - t, z - j19t))
            MPZ.mul!(tx, _bigbin(temp, j8 + t - j4, z - j147))
            MPZ.mul!(tx, _bigbin(temp, t + j4 - j8, z - j789))
            MPZ.sub!(Pt, tx, Pt)
        end
        MPZ.mul!(ABC, Pt)
        if isodd(xh + yh + zh)
            MPZ.neg!(ABC)
        end
        MPZ.add!(sum, ABC)
    end
    return
end

# D = \tan^2\beta = (m_1\omega_1)/(m_2\omega_2)
function _Moshinsky(N::Int, L::Int, n::Int, l::Int, n1::Int, l1::Int, n2::Int, l2::Int, Λ::Int, D::Rational{Int})
    check_Moshinsky(N, L, n, l, n1, l1, n2, l2, Λ) || return zero(SqrtRational{BigInt})
    e1 = 2 * n1 + l1
    e2 = 2 * n2 + l2
    E = 2 * N + L
    e = 2 * n + l
    χ = e1 + e2
    nl1 = n1 + l1
    nl2 = n2 + l2
    NL = N + L
    nl = n + l
    numD = numerator(D)
    denD = denominator(D)

    temp = BigInt()
    tx = BigInt()
    Pt = BigInt()
    ABC = BigInt()
    nu_a = BigInt()
    de_a = BigInt()
    nu_b = BigInt()
    de_b = BigInt()
    nu_c = BigInt()
    de_c = BigInt()
    nu_d = BigInt()
    de_d = BigInt()
    nu_fa = BigInt()
    de_fa = BigInt()
    M9j = BigInt()

    Snu = BigInt()
    Sde = BigInt()

    _bigbin(Snu, χ + 2, e1 + 1)
    _bigbin(Sde, χ + 2, E + 1)
    MPZ.mul!(Sde, _bigbin(temp, L + l + Λ + 1, 2Λ + 1))
    MPZ.mul!(Sde, _bigbin(temp, 2Λ, Λ + L - l))
    MPZ.mul!(Sde, _bigbin(temp, l1 + l2 + Λ + 1, 2Λ + 1))
    MPZ.mul!(Sde, _bigbin(temp, 2Λ, Λ + l1 - l2))
    _divgcd!(temp, Snu, Sde)

    MPZ.mul!(Snu, _bigbin(temp, 2nl1 + 1, nl1))
    MPZ.mul_ui!(Sde, convert(Culong, 2l1 + 1))
    MPZ.mul!(Sde, _bigbin(temp, e1 + 1, n1))

    MPZ.mul!(Snu, _bigbin(temp, 2nl2 + 1, nl2))
    MPZ.mul_ui!(Sde, convert(Culong, 2l2 + 1))
    MPZ.mul!(Sde, _bigbin(temp, e2 + 1, n2))

    MPZ.mul!(Snu, _bigbin(temp, 2NL + 1, NL))
    MPZ.mul_ui!(Sde, convert(Culong, 2L + 1))
    MPZ.mul!(Sde, _bigbin(temp, E + 1, N))

    MPZ.mul!(Snu, _bigbin(temp, 2nl + 1, nl))
    MPZ.mul_ui!(Sde, convert(Culong, 2l + 1))
    MPZ.mul!(Sde, _bigbin(temp, e + 1, n))

    _divgcd!(temp, Snu, Sde)

    MPZ.set_ui!(tx, convert(Culong, e1 + 2))
    MPZ.mul_ui!(tx, convert(Culong, e2 + 2))
    MPZ.mul_ui!(tx, convert(Culong, 2Λ + 1))
    MPZ.pow_ui!(tx, 2)
    MPZ.mul_2exp!(tx, l1 + l2 + L + l)
    MPZ.mul!(Sde, tx)
    _divgcd!(temp, Snu, Sde)

    if e - e1 > 0
        __gmpz_ui_pow_ui!(temp, convert(Culong, numD), convert(Culong, e - e1))
        MPZ.mul!(Snu, temp)
    elseif e - e1 < 0
        __gmpz_ui_pow_ui!(temp, convert(Culong, numD), convert(Culong, e1 - e))
        MPZ.mul!(Sde, temp)
    end
    __gmpz_ui_pow_ui!(temp, convert(Culong, denD), convert(Culong, e1 + E))
    MPZ.mul!(Snu, temp)
    __gmpz_ui_pow_ui!(temp, convert(Culong, numD + denD), convert(Culong, e + E))
    MPZ.mul!(Sde, temp)
    _divgcd!(temp, Snu, Sde)

    sum = zero(Rational{BigInt})
    for ea = 0:min(e1, E)
        eb = e1 - ea
        ec = E - ea
        ed = e2 - ec
        ed < 0 && continue
        _bigbin(nu_fa, e1 + 2, ea + 1)
        MPZ.mul!(nu_fa, _bigbin(temp, e2 + 2, ec + 1))
        __gmpz_ui_pow_ui!(temp, convert(Culong, numD), convert(Culong, ea))
        MPZ.mul!(nu_fa, temp)
        __gmpz_ui_pow_ui!(de_fa, convert(Culong, denD), convert(Culong, ea))
        _divgcd!(temp, nu_fa, de_fa)
        for la = rem(ea, 2):2:ea
            na = div(ea - la, 2)
            nla = na + la
            MPZ.mul_ui!(nu_a, nu_fa, convert(Culong, 2la + 1))
            MPZ.mul!(nu_a, _bigbin(temp, ea + 1, na))
            MPZ.mul_2exp!(nu_a, la)
            MPZ.mul!(de_a, de_fa, _bigbin(temp, 2nla + 1, nla))
            _divgcd!(temp, nu_a, de_a)
            for lb = abs(l1 - la):2:min(l1 + la, eb)
                nb = div(eb - lb, 2)
                nlb = nb + lb
                MPZ.mul_ui!(nu_b, nu_a, convert(Culong, 2lb + 1))
                MPZ.mul!(nu_b, _bigbin(temp, eb + 1, nb))
                MPZ.mul_2exp!(nu_b, lb)
                MPZ.mul!(de_b, de_a, _bigbin(temp, 2nlb + 1, nlb))
                MPZ.mul!(nu_b, _omega(tx, temp, la, l1, div(la + lb + l1, 2)))
                # Δ(lalbl1)
                MPZ.mul!(de_b, _bigbin(temp, la + lb + l1 + 1, 2l1 + 1))
                MPZ.mul!(de_b, _bigbin(temp, 2l1, l1 + la - lb))
                _divgcd!(temp, nu_b, de_b)
                for lc = abs(L - la):2:min(L + la, ec)
                    nc = div(ec - lc, 2)
                    nlc = nc + lc
                    MPZ.mul_ui!(nu_c, nu_b, convert(Culong, 2lc + 1))
                    MPZ.mul!(nu_c, _bigbin(temp, ec + 1, nc))
                    MPZ.mul_2exp!(nu_c, lc)
                    MPZ.mul!(de_c, de_b, _bigbin(temp, 2nlc + 1, nlc))
                    MPZ.mul!(nu_c, _omega(tx, temp, la, L, div(la + lc + L, 2)))
                    # Δ(lalcl3)
                    MPZ.mul!(de_c, _bigbin(temp, la + lc + L + 1, 2L + 1))
                    MPZ.mul!(de_c, _bigbin(temp, 2L, L + la - lc))
                    _divgcd!(temp, nu_c, de_c)
                    ldmin = max(abs(l2 - lc), abs(l - lb))
                    ldmax = min(ed, l2 + lc, l + lb)
                    for ld = ldmin:2:ldmax
                        nd = div(ed - ld, 2)
                        nld = nd + ld
                        MPZ.mul_ui!(nu_d, nu_c, convert(Culong, 2ld + 1))
                        MPZ.mul!(nu_d, _bigbin(temp, ed + 1, nd))
                        MPZ.mul_2exp!(nu_d, ld)
                        MPZ.mul!(de_d, de_c, _bigbin(temp, 2nld + 1, nld))
                        MPZ.mul!(nu_d, _omega(tx, temp, lc, l2, div(lc + ld + l2, 2)))
                        MPZ.mul!(nu_d, _omega(tx, temp, lb, l, div(lb + ld + l, 2)))
                        # Δ(lbldl4)
                        MPZ.mul!(de_d, _bigbin(temp, lb + ld + l + 1, 2l + 1))
                        MPZ.mul!(de_d, _bigbin(temp, 2l, l + lb - ld))
                        # Δ(lcldl2)
                        MPZ.mul!(de_d, _bigbin(temp, lc + ld + l2 + 1, 2l2 + 1))
                        MPZ.mul!(de_d, _bigbin(temp, 2l2, l2 + lc - ld))
                        _divgcd!(temp, nu_d, de_d)
                        if isodd(ld)
                            MPZ.neg!(nu_d)
                        end
                        _m9j!(M9j, temp, tx, Pt, ABC, la, lb, l1, lc, ld, l2, L, l, Λ)
                        # don't use MPQ.add! to avoid memory allocation
                        MPZ.mul!(nu_d, M9j)
                        _divgcd!(temp, nu_d, de_d)
                        __gmpq_add!(tx, sum, nu_d, de_d)
                    end
                end
            end
        end
    end
    simplify4!(tx, sum.num, sum.den, Snu, Sde, max(χ + 2, max(l1 + l2, L + l) + Λ + 1, numD + denD))
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
    xGaunt(l1::Integer, l2::Integer, l3::Integer, m1::Integer, m2::Integer, m3::Integer)
Gaunt coefficient with out ``1/\sqrt{\pi}`` factor.
(This package can only handle integer arithmetic, it is easy to add the factor mannually.)
```math
\begin{aligned}
\mathrm{Gaunt}(l_1 l_2 l_3 m_1 m_2 m_3) &= \int Y_{l_1 m_1}(\theta, \phi) Y_{l_2 m_2}(\theta, \phi) Y_{l_3 m_3}(\theta, \phi) d\Omega \\
&= \sqrt{\frac{(2l_1 + 1)(2l_2 + 1)(2l_3 + 1)}{4\pi}} \begin{pmatrix} l_1 & l_2 & l_3 \\ m_1 & m_2 & m_3 \end{pmatrix} \begin{pmatrix} l_1 & l_2 & l_3 \\ 0 & 0 & 0 \end{pmatrix}. \\
\end{aligned}
```
"""
@inline xGaunt(l1::Integer, l2::Integer, l3::Integer, m1::Integer, m2::Integer, m3::Integer) = _xGaunt(Int.((l1, l2, l3, m1, m2, m3))...)

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
    Moshinsky(N::Integer, L::Integer, n::Integer, l::Integer, n1::Integer, l1::Integer, n2::Integer, l2::Integer, Λ::Integer, D::Union{Rational,Integer}=1)
Moshinsky bracket，Ref: Buck et al. Nuc. Phys. A 600 (1996) 387-402.
This function is designed for demonstration the exact result,
It only calculate the case of ``\tan(\beta) = \sqrt{m_1\omega_1/m_2\omega_2} = 1``.
"""
@inline Moshinsky(N::Integer, L::Integer, n::Integer, l::Integer, n1::Integer, l1::Integer, n2::Integer, l2::Integer, Λ::Integer, D::Union{Rational,Integer}=1) = begin
    iN, iL, in, il, in1, il1, in2, il2, iΛ = Int.((N, L, n, l, n1, l1, n2, l2, Λ))
    numD = Int(numerator(D))
    denD = Int(denominator(D))
    return _Moshinsky(iN, iL, in, il, in1, il1, in2, il2, iΛ, Rational(numD, denD))
end
