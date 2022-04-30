# This file contains the core functions of WignerSymbols and CG-coefficient.

"""
    HalfInt = Union{Integer, Rational}
Angular momentum quantum numbers may be half integers and integers. With `HalfInt` type, you can use both `2` and `3//2` as parameters.
But the parameter like `5//3`, who's denominator is not `2` while gives out error.
"""
const HalfInt = Union{Integer,Rational}

# basic CG coefficient calculation function
function _dCG(dj1::Int64, dj2::Int64, dj3::Int64, dm1::Int64, dm2::Int64, dm3::Int64)
    check_jm(dj1, dm1) || error(jm_mismatching_msg(dj1, dm1))
    check_jm(dj2, dm2) || error(jm_mismatching_msg(dj2, dm2))
    check_jm(dj3, dm3) || error(jm_mismatching_msg(dj3, dm3))
    check_couple(dj1, dj2, dj3) || error(miss_couple_msg(dj1, dj2, dj3))
    dm1 + dm2 == dm3 || return zero(SqrtRational)
    J::Int64 = div(dj1 + dj2 + dj3, 2)
    Jm1::Int64 = J - dj1
    Jm2::Int64 = J - dj2
    Jm3::Int64 = J - dj3
    j1mm1::Int64 = div(dj1 - dm1, 2)
    j2mm2::Int64 = div(dj2 - dm2, 2)
    j3mm3::Int64 = div(dj3 - dm3, 2)
    j2pm2::Int64 = div(dj2 + dm2, 2)
    # value in sqrt
    A = widemul(binomial(dj1, Jm2), binomial(dj2, Jm3)) // widemul(
            widemul(binomial(J + 1, Jm3), binomial(dj1, j1mm1)),
            widemul(binomial(dj2, j2mm2), binomial(dj3, j3mm3))
        )
    B::BigInt = zero(BigInt)
    low::Int64 = max(zero(Int64), j1mm1 - Jm2, j2pm2 - Jm1)
    high::Int64 = min(Jm3, j1mm1, j2pm2)
    for z in low:high
        B = -B + widemul(widemul(binomial(Jm3, z), binomial(Jm2, j1mm1 - z)), binomial(Jm1, j2pm2 - z))
    end
    return SqrtRational(iphase(high) * B, A)
end

# use CG coefficient to calculate 3j-symbol
function _d3j(dj1::Int64, dj2::Int64, dj3::Int64, dm1::Int64, dm2::Int64, dm3::Int64)
    iphase(dj1 + div(dj3 + dm3, 2)) * exact_sqrt(1 // (dj3 + 1)) * _dCG(dj1, dj2, dj3, -dm1, -dm2, dm3)
end

# basic 6j-symbol calculation funciton
function _d6j(dj1::Int64, dj2::Int64, dj3::Int64, dj4::Int64, dj5::Int64, dj6::Int64)
    check_couple(dj1, dj2, dj3) || error(miss_couple_msg(dj1, dj2, dj3))
    check_couple(dj1, dj5, dj6) || error(miss_couple_msg(dj1, dj5, dj6))
    check_couple(dj4, dj2, dj6) || error(miss_couple_msg(dj4, dj2, dj6))
    check_couple(dj4, dj5, dj3) || error(miss_couple_msg(dj4, dj5, dj3))
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
    # value in sqrt
    A = widemul(binomial(j123 + 1, dj1 + 1), binomial(dj1, jpm123)) // widemul(
            widemul(
                widemul(binomial(j156 + 1, dj1 + 1), binomial(dj1, jpm156)),
                widemul(binomial(j453 + 1, dj4 + 1), binomial(dj4, jpm453))
            ),
            widemul(binomial(j426 + 1, dj4 + 1), binomial(dj4, jpm426))
        )
    B::BigInt = zero(BigInt)
    low::Int64 = max(j123, j453, j426, j156)
    high::Int64 = min(jpm123 + j453, jpm132 + j426, jpm231 + j156)
    for x = low:high
        B = -B + widemul(
            widemul(binomial(x + 1, j123 + 1), binomial(jpm123, x - j453)),
            widemul(binomial(jpm132, x - j426), binomial(jpm231, x - j156))
        )
    end
    return SqrtRational(iphase(high) * B//(dj4+1), A)
end

# use 6j-symbol to calculate Racah coefficient
function _dRacah(dj1::Int64, dj2::Int64, dj3::Int64, dj4::Int64, dj5::Int64, dj6::Int64)
    iphase(div(dj1+dj2+dj3+dj4, 2)) * _d6j(dj1, dj2, dj5, dj4, dj3, dj6)
end

# basic 9j-symbol calculation funciton
function _d9j(dj1::Int64, dj2::Int64, dj3::Int64,
              dj4::Int64, dj5::Int64, dj6::Int64,
              dj7::Int64, dj8::Int64, dj9::Int64)
    check_couple(dj1, dj2, dj3) || error(miss_couple_msg(dj1, dj2, dj3))
    check_couple(dj4, dj5, dj6) || error(miss_couple_msg(dj4, dj5, dj6))
    check_couple(dj7, dj8, dj9) || error(miss_couple_msg(dj7, dj8, dj9))
    check_couple(dj1, dj4, dj7) || error(miss_couple_msg(dj1, dj4, dj7))
    check_couple(dj2, dj5, dj8) || error(miss_couple_msg(dj2, dj5, dj8))
    check_couple(dj3, dj6, dj9) || error(miss_couple_msg(dj3, dj6, dj9))
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
    P0_nu = widemul(widemul(
                widemul(binomial(j123 + 1, dj1 + 1), binomial(dj1, pm123)),
                widemul(binomial(j456 + 1, dj5 + 1), binomial(dj5, pm456))
                ),
                widemul(binomial(j789 + 1, dj9 + 1), binomial(dj9, pm798))
            )
    P0_de = widemul(widemul(
                widemul(binomial(j147 + 1, dj1 + 1), binomial(dj1, div(dj1 + dj4 - dj7, 2))),
                widemul(binomial(j258 + 1, dj5 + 1), binomial(dj5, div(dj2 + dj5 - dj8, 2)))
                ),
                widemul(binomial(j369 + 1, dj9 + 1), binomial(dj9, div(dj3 + dj9 - dj6, 2)))
            )
    P0 = P0_nu // P0_de
    dtl::Int64 = max(abs(dj2 - dj6), abs(dj4 - dj8), abs(dj1 - dj9))
    dth::Int64 = min(dj2 + dj6, dj4 + dj8, dj1 + dj9)
    PABC::Rational{BigInt} = zero(Rational{BigInt})
    for dt::Int64 = dtl:2:dth
        j19t::Int64 = div(dj1 + dj9 + dt, 2)
        j26t::Int64 = div(dj2 + dj6 + dt, 2)
        j48t::Int64 = div(dj4 + dj8 + dt, 2)
        Pt_de = widemul(widemul(
                    widemul(binomial(j19t + 1, dt + 1), binomial(dt, div(dj1 + dt - dj9, 2))),
                    widemul(binomial(j26t + 1, dt + 1), binomial(dt, div(dj2 + dt - dj6, 2)))
                    ),
                    widemul(binomial(j48t + 1, dt + 1), binomial(dt, div(dj4 + dt - dj8, 2)))
                )
        Pt_de *= (dt + 1)^2
        xl::Int64 = max(j123, j369, j26t, j19t)
        xh::Int64 = min(pm123 + j369, pm132 + j26t, pm231 + j19t)
        At::BigInt = zero(BigInt)
        for x = xl:xh
            At = -At + widemul(
                    widemul(binomial(x + 1, j123 + 1), binomial(pm123, x - j369)),
                    widemul(binomial(pm132, x - j26t), binomial(pm231, x - j19t))
                )
        end
        yl::Int64 = max(j456, j26t, j258, j48t)
        yh::Int64 = min(pm456 + j26t, pm465 + j258, pm564 + j48t)
        Bt::BigInt = zero(BigInt)
        for y = yl:yh
            Bt = -Bt + widemul(
                    widemul(binomial(y + 1, j456 + 1), binomial(pm456, y - j26t)),
                    widemul(binomial(pm465, y - j258), binomial(pm564, y - j48t))
                )
        end
        zl::Int64 = max(j789, j19t, j48t, j147)
        zh::Int64 = min(pm789 + j19t, pm798 + j48t, pm897 + j147)
        Ct::BigInt = zero(BigInt)
        for z = zl:zh
            Ct = -Ct + widemul(
                    widemul(binomial(z + 1, j789 + 1), binomial(pm789, z - j19t)),
                    widemul(binomial(pm798, z - j48t), binomial(pm897, z - j147))
                )
        end
        PABC += (iphase(xh + yh + zh) * At * Bt * Ct) // Pt_de
    end
    return SqrtRational(iphase(dth) * PABC, P0)
end

"""
    dCG(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer)
CG coefficient function with double angular monentum number parameters, so that the parameters can be integer.
You can calculate `dCG(1, 1, 2, 1, 1, 2)` to calculate the real `CG(1//2, 1//2, 1, 1/2, 1//2, 1)`
"""
dCG(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer) = _dCG(Int64.((dj1, dj2, dj3, dm1, dm2, dm3))...)

"""
    d3j(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer)
3j-symbol function with double angular monentum number parameters, so that the parameters can be integer.
"""
d3j(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer) = _d3j(Int64.((dj1, dj2, dj3, dm1, dm2, dm3))...)

"""
    d6j(dj1::Integer, dj2::Integer, dj3::Integer, dj4::Integer, dj5::Integer, dj6::Integer)
6j-symbol function with double angular monentum number parameters, so that the parameters can be integer.
"""
d6j(dj1::Integer, dj2::Integer, dj3::Integer, dj4::Integer, dj5::Integer, dj6::Integer) = _d6j(Int64.((dj1, dj2, dj3, dj4, dj5, dj6))...)

"""
    dRacah(dj1::Integer, dj2::Integer, dj3::Integer, dj4::Integer, dj5::Integer, dj6::Integer)
Racah coefficient function with double angular momentum parameters, so that the parameters can be integer.
"""
dRacah(dj1::Integer, dj2::Integer, dj3::Integer, dj4::Integer, dj5::Integer, dj6::Integer) = _dRacah(Int64.((dj1, dj2, dj3, dj4, dj5, dj6))...)

"""
    d9j(dj1::Integer, dj2::Integer, dj3::Integer,
        dj4::Integer, dj5::Integer, dj6::Integer,
        dj7::Integer, dj8::Integer, dj9::Integer)
9j-symbol function with double angular monentum number parameters, so that the parameters can be integer.
"""
d9j(dj1::Integer, dj2::Integer, dj3::Integer,
    dj4::Integer, dj5::Integer, dj6::Integer,
    dj7::Integer, dj8::Integer, dj9::Integer) = _d9j(Int64.((dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9))...)

@doc raw"""
    CG(j1::HalfInt, j2::HalfInt, j3::HalfInt, m1::HalfInt, m2::HalfInt, m3::HalfInt)
CG coefficient ``\langle j_1m_1 j_2m_2 | j_3m_3 \rangle``
"""
CG(j1::HalfInt, j2::HalfInt, j3::HalfInt, m1::HalfInt, m2::HalfInt, m3::HalfInt) = _dCG(Int64.((2j1, 2j2, 2j3, 2m1, 2m2, 2m3))...)

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
threeJ(j1::HalfInt, j2::HalfInt, j3::HalfInt, m1::HalfInt, m2::HalfInt, m3::HalfInt) = _d3j(Int64.((2j1, 2j2, 2j3, 2m1, 2m2, 2m3))...)

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
sixJ(j1::HalfInt, j2::HalfInt, j3::HalfInt, j4::HalfInt, j5::HalfInt, j6::HalfInt) = _d6j(Int64.((2j1, 2j2, 2j3, 2j4, 2j5, 2j6))...)

@doc raw"""
    Racah(j1::HalfInt, j2::HalfInt, j3::HalfInt, j4::HalfInt, j5::HalfInt, j6::HalfInt)
Racah coefficient
```math
W(j_1j_2j_3j_4, j_5j_6) = (-1)^{j_1+j_2+j_3+j_4} \begin{Bmatrix}
j_1 & j_2 & j_5 \\
j_4 & j_3 & j_6
\end{Bmatrix}
```
"""
Racah(j1::HalfInt, j2::HalfInt, j3::HalfInt, j4::HalfInt, j5::HalfInt, j6::HalfInt) = _dRacah(Int64.((2j1, 2j2, 2j3, 2j4, 2j5, 2j6))...)

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
      j7::HalfInt, j8::HalfInt, j9::HalfInt) = _d9j(Int64.((2j1, 2j2, 2j3, 2j4, 2j5, 2j6, 2j7, 2j8, 2j9))...)
