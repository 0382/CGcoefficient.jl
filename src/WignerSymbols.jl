# This file contains the core functions of WignerSymbols and CG-coefficient. 

"""
    HalfInt = Union{Integer, Rational}
Angular momentum quantum numbers may be half integers and integers. With `HalfInt` type, you can use both `2` and `3//2` as parameters.
But the parameter like `5//3`, who's denominator is not `2` while gives out error.
"""
const HalfInt = Union{Integer,Rational}

"""
basic CG coefficient calculation function
"""
function _dCG(dj1::BigInt, dj2::BigInt, dj3::BigInt, dm1::BigInt, dm2::BigInt, dm3::BigInt)
    check_jm(dj1, dm1) || error(jm_mismatching_msg(dj1, dm1))
    check_jm(dj2, dm2) || error(jm_mismatching_msg(dj2, dm2))
    check_jm(dj3, dm3) || error(jm_mismatching_msg(dj3, dm3))
    check_couple(dj1, dj2, dj3) || error(miss_couple_msg(dj1, dj2, dj3))
    dm1 + dm2 == dm3 || return zero(SqrtRational)
    J::BigInt = div(dj1 + dj2 + dj3, 2)
    Jm1::BigInt = J - dj1
    Jm2::BigInt = J - dj2
    Jm3::BigInt = J - dj3
    j1mm1::BigInt = div(dj1 - dm1, 2)
    j2mm2::BigInt = div(dj2 - dm2, 2)
    j3mm3::BigInt = div(dj3 - dm3, 2)
    j2pm2::BigInt = div(dj2 + dm2, 2)
    A = SqrtRational(
        binomial(dj1, Jm2) * binomial(dj2, Jm3) // (
            binomial(J + 1, Jm3) * binomial(dj1, j1mm1) * binomial(dj2, j2mm2) * binomial(dj3, j3mm3)
        )
    )
    B::BigInt = zero(BigInt)
    low::BigInt = max(zero(BigInt), j1mm1 - Jm2, j2pm2 - Jm1)
    high::BigInt = min(Jm3, j1mm1, j2pm2)
    for z in low:high
        B = -B + binomial(Jm3, z) * binomial(Jm2, j1mm1 - z) * binomial(Jm1, j2pm2 - z)
    end
    return A * (iphase(high) * B)
end

"""
use CG coefficient to calculate 3j-symbol
"""
function _d3j(dj1::BigInt, dj2::BigInt, dj3::BigInt, dm1::BigInt, dm2::BigInt, dm3::BigInt)
    iphase(dj1 + div(dj3 + dm3, 2)) * SqrtRational(1 // (dj3 + 1)) * _dCG(dj1, dj2, dj3, -dm1, -dm2, dm3)
end

"""
basic 6j-symbol calculation funciton
"""
function _d6j(dj1::BigInt, dj2::BigInt, dj3::BigInt, dj4::BigInt, dj5::BigInt, dj6::BigInt)
    check_couple(dj1, dj2, dj3) || error(miss_couple_msg(dj1, dj2, dj3))
    check_couple(dj1, dj5, dj6) || error(miss_couple_msg(dj1, dj5, dj6))
    check_couple(dj4, dj2, dj6) || error(miss_couple_msg(dj4, dj2, dj6))
    check_couple(dj4, dj5, dj3) || error(miss_couple_msg(dj4, dj5, dj3))
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
    A = SqrtRational(
        binomial(j123 + 1, dj1 + 1) * binomial(dj1, jpm123) // (
            binomial(j156 + 1, dj1 + 1) * binomial(dj1, jpm156) *
            binomial(j453 + 1, dj4 + 1) * binomial(dj4, jpm453) *
            binomial(j426 + 1, dj4 + 1) * binomial(dj4, jpm426)
        )
    )
    A = A * Rational{BigInt}(1 // (dj4 + 1))
    B::BigInt = zero(BigInt)
    low::BigInt = max(j123, j453, j426, j156)
    high::BigInt = min(jpm123 + j453, jpm132 + j426, jpm231 + j156)
    for n = low:high
        B = -B + binomial(n + 1, j123 + 1) * binomial(jpm123, n - j453) *
                binomial(jpm132, n - j426) * binomial(jpm231, n - j156)
    end
    return A * (iphase(high) * B)
end

"""
basic 9j-symbol calculation funciton
"""
function _d9j(dj1::BigInt, dj2::BigInt, dj3::BigInt,
              dj4::BigInt, dj5::BigInt, dj6::BigInt,
              dj7::BigInt, dj8::BigInt, dj9::BigInt)
    check_couple(dj1, dj2, dj3) || error(miss_couple_msg(dj1, dj2, dj3))
    check_couple(dj4, dj5, dj6) || error(miss_couple_msg(dj4, dj5, dj6))
    check_couple(dj7, dj8, dj9) || error(miss_couple_msg(dj7, dj8, dj9))
    check_couple(dj1, dj4, dj7) || error(miss_couple_msg(dj1, dj4, dj7))
    check_couple(dj2, dj5, dj8) || error(miss_couple_msg(dj2, dj5, dj8))
    check_couple(dj3, dj6, dj9) || error(miss_couple_msg(dj3, dj6, dj9))
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
    P0_nu::BigInt = binomial(j123 + 1, dj1 + 1) * binomial(dj1, pm123) *
                    binomial(j456 + 1, dj5 + 1) * binomial(dj5, pm456) *
                    binomial(j789 + 1, dj9 + 1) * binomial(dj9, pm798)
    P0_de::BigInt = binomial(j147 + 1, dj1 + 1) * binomial(dj1, div(dj1 + dj4 - dj7, 2)) *
                    binomial(j258 + 1, dj5 + 1) * binomial(dj5, div(dj2 + dj5 - dj8, 2)) *
                    binomial(j369 + 1, dj9 + 1) * binomial(dj9, div(dj3 + dj9 - dj6, 2))
    P0 = SqrtRational(P0_nu // P0_de)
    dtl::BigInt = max(abs(dj2 - dj6), abs(dj4 - dj8), abs(dj1 - dj9))
    dth::BigInt = min(dj2 + dj6, dj4 + dj8, dj1 + dj9)
    PABC::Rational{BigInt} = zero(Rational{BigInt})
    for dt::BigInt = dtl:2:dth
        j19t::BigInt = div(dj1 + dj9 + dt, 2)
        j26t::BigInt = div(dj2 + dj6 + dt, 2)
        j48t::BigInt = div(dj4 + dj8 + dt, 2)
        Pt_de::BigInt = binomial(j19t + 1, dt + 1) * binomial(dt, div(dj1 + dt - dj9, 2)) *
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
    return P0 * (iphase(dth) * PABC)
end

"""
    dCG(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer)
CG coefficient function with double angular monentum number parameters, so that the parameters can be integer.
You can calculate `dCG(1, 1, 2, 1, 1, 2)` to calculate the real `CG(1//2, 1//2, 1, 1/2, 1//2, 1)`
"""
dCG(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer) = _dCG(BigInt.((dj1, dj2, dj3, dm1, dm2, dm3))...)

"""
    d3j(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer)
3j-symbol function with double angular monentum number parameters, so that the parameters can be integer.
"""
d3j(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer) = _d3j(BigInt.((dj1, dj2, dj3, dm1, dm2, dm3))...)

"""
    d6j(dj1::Integer, dj2::Integer, dj3::Integer, dj4::Integer, dj5::Integer, dj6::Integer)
6j-symbol function with double angular monentum number parameters, so that the parameters can be integer.
"""
d6j(dj1::Integer, dj2::Integer, dj3::Integer, dj4::Integer, dj5::Integer, dj6::Integer) = _d6j(BigInt.((dj1, dj2, dj3, dj4, dj5, dj6))...)

"""
    d9j(dj1::Integer, dj2::Integer, dj3::Integer,
        dj4::Integer, dj5::Integer, dj6::Integer,
        dj7::Integer, dj8::Integer, dj9::Integer)
9j-symbol function with double angular monentum number parameters, so that the parameters can be integer.
"""
d9j(dj1::Integer, dj2::Integer, dj3::Integer,
    dj4::Integer, dj5::Integer, dj6::Integer,
    dj7::Integer, dj8::Integer, dj9::Integer) = _d9j(BigInt.((dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9))...)

@doc raw"""
    CG(j1::HalfInt, j2::HalfInt, j3::HalfInt, m1::HalfInt, m2::HalfInt, m3::HalfInt)
CG coefficient ``\langle j_1m_1 j_2m_2 | j_3m_3 \rangle``
"""
CG(j1::HalfInt, j2::HalfInt, j3::HalfInt, m1::HalfInt, m2::HalfInt, m3::HalfInt) = _dCG(BigInt.((2j1, 2j2, 2j3, 2m1, 2m2, 2m3))...)

@doc raw"""
    threeJ(j1::HalfInt, j2::HalfInt, j3::HalfInt, m1::HalfInt, m2::HalfInt, m3::HalfInt)
Wigner 3j-symbol
```math
\begin{pmatrix}
j_1 & j_2 & j_3 \\
m_1 & m_2 & m_3
\end{pmatrix}
```
"""
threeJ(j1::HalfInt, j2::HalfInt, j3::HalfInt, m1::HalfInt, m2::HalfInt, m3::HalfInt) = _d3j(BigInt.((2j1, 2j2, 2j3, 2m1, 2m2, 2m3))...)

@doc raw"""
    sixJ(j1::HalfInt, j2::HalfInt, j3::HalfInt, j4::HalfInt, j5::HalfInt, j6::HalfInt)
Wigner 6j-symbol
```math
\begin{Bmatrix}
j_1 & j_2 & j_3 \\
j_4 & j_5 & j_6
\end{Bmatrix}
```
"""
sixJ(j1::HalfInt, j2::HalfInt, j3::HalfInt, j4::HalfInt, j5::HalfInt, j6::HalfInt) = _d6j(BigInt.((2j1, 2j2, 2j3, 2j4, 2j5, 2j6))...)

@doc raw"""
    nineJ(j1::HalfInt, j2::HalfInt, j3::HalfInt,
          j4::HalfInt, j5::HalfInt, j6::HalfInt,
          j7::HalfInt, j8::HalfInt, j9::HalfInt)
Wigner 9j-symbol
```math
\begin{Bmatrix}
j_1 & j_2 & j_3 \\
j_4 & j_5 & j_6 \\
j_7 & j_8 & j_9
\end{Bmatrix}
```
"""
nineJ(j1::HalfInt, j2::HalfInt, j3::HalfInt,
      j4::HalfInt, j5::HalfInt, j6::HalfInt,
      j7::HalfInt, j8::HalfInt, j9::HalfInt) = _d9j(BigInt.((2j1, 2j2, 2j3, 2j4, 2j5, 2j6, 2j7, 2j8, 2j9))...)
