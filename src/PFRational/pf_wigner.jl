function stagger_sum!(sum::BigInt, cf::PFRational{T}, xs::Vector{PFRational{T}}) where T<:Integer
    Base.GMP.MPZ.set_ui!(sum, 0)
    _copy!(cf, xs[1])
    for i in 2:length(xs)
        _gcd!(cf, xs[i])
    end
    t = zero(BigInt)
    for i in 1:length(xs)
        _div!(xs[i], cf)
        _numerator!(t, xs[i])
        if isodd(i)
            Base.GMP.MPZ.add!(sum, t)
        else
            Base.GMP.MPZ.sub!(sum, t)
        end
    end
end

const c_tdiv_q_ui = Base.GMP.MPZ.gmpz(:tdiv_q_ui)

@inline function __gmpz_tdiv_q_ui!(q::BigInt, n::BigInt, d::Culong)
    ccall(c_tdiv_q_ui, Culong, (Ref{BigInt}, Ref{BigInt}, Culong), q, n, d)
end

function extract_to!(pf::PFRational{T}, num::BigInt) where T<:Integer
    (iszero(num) || isone(abs(num))) && return
    t = zero(BigInt)
    for i in eachindex(pf.e)
        r = __gmpz_tdiv_q_ui!(t, num, nth_stored_prime(i))
        while iszero(r)
            pf.e[i] += 1
            Base.GMP.MPZ.set!(num, t)
            r = __gmpz_tdiv_q_ui!(t, num, nth_stored_prime(i))
        end
    end
    return
end

# assume `sd == rn == rd == 1`
# s * √r
function extract_sr!(sn::BigInt, sd::BigInt, rn::BigInt, rd::BigInt, pf::PFRational{T}) where T<:Integer
    for i in eachindex(pf.e)
        si, ri = divrem(pf.e[i], 2)
        if si > 0
            mul_p_pow_k!(sn, i, si)
        elseif si < 0
            mul_p_pow_k!(sd, i, -si)
        end
        if ri > 0
            Base.GMP.MPZ.mul_ui!(rn, nth_stored_prime(i))
        elseif ri < 0
            Base.GMP.MPZ.mul_ui!(rd, nth_stored_prime(i))
        end
    end
    return
end

# B * √A
function _pf_CG_impl!(A::PFRational{Int16}, B::BigInt, dj1::Int, dj2::Int, dj3::Int, dm1::Int, dm2::Int, dm3::Int)
    J::Int = div(dj1 + dj2 + dj3, 2)
    Jm1::Int = J - dj1
    Jm2::Int = J - dj2
    Jm3::Int = J - dj3
    j1mm1::Int = div(dj1 - dm1, 2)
    j2mm2::Int = div(dj2 - dm2, 2)
    j3mm3::Int = div(dj3 - dm3, 2)
    j2pm2::Int = div(dj2 + dm2, 2)

    Jsize::Int = max_prime_index(J + 1)
    low::Int = max(zero(Int), j1mm1 - Jm2, j2pm2 - Jm1)
    high::Int = min(Jm3, j1mm1, j2pm2)

    pf_alloc!(A, Jsize)
    if low == high
        z = low
        _copy!(A, _pf_bin(Jm3, z))
        _mul!(A, _pf_bin(Jm2, j1mm1 - z))
        _mul!(A, _pf_bin(Jm1, j2pm2 - z))
        Base.GMP.MPZ.set_ui!(B, 1)
    else
        Bs = Vector{PFRational{Int16}}(undef, high - low + 1)
        for z in low:high
            Bs[z - low + 1] = pf_alloc(Int16, Jsize)
            _copy!(Bs[z - low + 1], _pf_bin(Jm3, z))
            _mul!(Bs[z - low + 1], _pf_bin(Jm2, j1mm1 - z))
            _mul!(Bs[z - low + 1], _pf_bin(Jm1, j2pm2 - z))
        end
        stagger_sum!(B, A, Bs)
        iszero(B) && return
        extract_to!(A, B)
    end
    square!(A)
    _mul!(A, _pf_bin(dj1, Jm2))
    _mul!(A, _pf_bin(dj2, Jm3))
    _div!(A, _pf_bin(J + 1, Jm3))
    _div!(A, _pf_bin(dj1, j1mm1))
    _div!(A, _pf_bin(dj2, j2mm2))
    _div!(A, _pf_bin(dj3, j3mm3))
    if isodd(low)
        Base.GMP.MPZ.neg!(B)
    end
    return
end

function _pf_3j_impl!(A::PFRational{Int16}, B::BigInt, dj1::Int, dj2::Int, dj3::Int, dm1::Int, dm2::Int, dm3::Int)
    J::Int = div(dj1 + dj2 + dj3, 2)
    Jm1::Int = J - dj1
    Jm2::Int = J - dj2
    Jm3::Int = J - dj3
    j1mm1::Int = div(dj1 - dm1, 2)
    j2mm2::Int = div(dj2 - dm2, 2)
    j3mm3::Int = div(dj3 - dm3, 2)
    j1pm1::Int = div(dj1 + dm1, 2)

    Jsize::Int = max_prime_index(J + 1)
    low::Int = max(zero(Int), j1pm1 - Jm2, j2mm2 - Jm1)
    high::Int = min(Jm3, j1pm1, j2mm2)
    
    pf_alloc!(A, Jsize)
    if low == high
        z = low
        _copy!(A, _pf_bin(Jm3, z))
        _mul!(A, _pf_bin(Jm2, j1pm1 - z))
        _mul!(A, _pf_bin(Jm1, j2mm2 - z))
        Base.GMP.MPZ.set_ui!(B, 1)
    else
        Bs = Vector{PFRational{Int16}}(undef, high - low + 1)
        for z in low:high
            Bs[z - low + 1] = pf_alloc(Int16, Jsize)
            _copy!(Bs[z - low + 1], _pf_bin(Jm3, z))
            _mul!(Bs[z - low + 1], _pf_bin(Jm2, j1pm1 - z))
            _mul!(Bs[z - low + 1], _pf_bin(Jm1, j2mm2 - z))
        end
        stagger_sum!(B, A, Bs)
        iszero(B) && return
        extract_to!(A, B)
    end
    square!(A)
    _mul!(A, _pf_bin(dj1, Jm2))
    _mul!(A, _pf_bin(dj2, Jm1))
    _div!(A, pf_int(J + 1))
    _div!(A, _pf_bin(J, Jm3))
    _div!(A, _pf_bin(dj1, j1mm1))
    _div!(A, _pf_bin(dj2, j2mm2))
    _div!(A, _pf_bin(dj3, j3mm3))
    if isodd(low + dj1 + div(dj3 + dm3, 2))
        Base.GMP.MPZ.neg!(B)
    end
    return
end

function _pf_6j_impl!(A::PFRational{Int16}, B::BigInt, dj1::Int, dj2::Int, dj3::Int, dj4::Int, dj5::Int, dj6::Int)
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
    low::Int = max(j123, j453, j426, j156)
    high::Int = min(jpm123 + j453, jpm132 + j426, jpm231 + j156)
    Jsize::Int = max_prime_index(high + 1)
    pf_alloc!(A, Jsize)
    if low == high
        x = low
        _copy!(A, _pf_bin(x + 1, j123 + 1))
        _mul!(A, _pf_bin(jpm123, x - j453))
        _mul!(A, _pf_bin(jpm132, x - j426))
        _mul!(A, _pf_bin(jpm231, x - j156))
        Base.GMP.MPZ.set_ui!(B, 1)
    else
        Bs = Vector{PFRational{Int16}}(undef, high - low + 1)
        for x in low:high
            Bs[x - low + 1] = pf_alloc(Int16, Jsize)
            _copy!(Bs[x - low + 1], _pf_bin(x + 1, j123 + 1))
            _mul!(Bs[x - low + 1], _pf_bin(jpm123, x - j453))
            _mul!(Bs[x - low + 1], _pf_bin(jpm132, x - j426))
            _mul!(Bs[x - low + 1], _pf_bin(jpm231, x - j156))
        end
        stagger_sum!(B, A, Bs)
        iszero(B) && return
        extract_to!(A, B)
    end
    _div!(A, pf_int(dj4 + 1))
    square!(A)
    _mul!(A, _pf_bin(j123 + 1, dj1 + 1))
    _mul!(A, _pf_bin(dj1, jpm123))
    _div!(A, _pf_bin(j156 + 1, dj1 + 1))
    _div!(A, _pf_bin(dj1, jpm156))
    _div!(A, _pf_bin(j453 + 1, dj4 + 1))
    _div!(A, _pf_bin(dj4, jpm453))
    _div!(A, _pf_bin(j426 + 1, dj4 + 1))
    _div!(A, _pf_bin(dj4, jpm426))
    if isodd(low)
        Base.GMP.MPZ.neg!(B)
    end
    return
end

function _pf_CG(dj1::Int, dj2::Int, dj3::Int, dm1::Int, dm2::Int, dm3::Int)
    check_CG(dj1, dj2, dj3, dm1, dm2, dm3) || return zero(SqrtRational{BigInt})
    A = one(PFRational{Int16})
    sn = one(BigInt)
    _pf_CG_impl!(A, sn, dj1, dj2, dj3, dm1, dm2, dm3)
    iszero(sn) && return zero(SqrtRational{BigInt})
    sd = one(BigInt)
    rn = one(BigInt)
    rd = one(BigInt)
    extract_sr!(sn, sd, rn, rd, A)
    return SqrtRational(Base.unsafe_rational(sn, sd), Base.unsafe_rational(rn, rd))
end

function _pf_fCG(dj1::Int, dj2::Int, dj3::Int, dm1::Int, dm2::Int, dm3::Int)::Float64
    check_CG(dj1, dj2, dj3, dm1, dm2, dm3) || return zero(Float64)
    A = one(PFRational{Int16})
    sn = one(BigInt)
    _pf_CG_impl!(A, sn, dj1, dj2, dj3, dm1, dm2, dm3)
    iszero(sn) && return zero(Float64)
    sd = one(BigInt)
    rn = one(BigInt)
    rd = one(BigInt)
    extract_sr!(sn, sd, rn, rd, A)
    ans = BigFloat(sn) / BigFloat(sd) * sqrt(BigFloat(rn) / BigFloat(rd))
    return convert(Float64, ans)
end

function _pf_3j(dj1::Int, dj2::Int, dj3::Int, dm1::Int, dm2::Int, dm3::Int)
    check_3j(dj1, dj2, dj3, dm1, dm2, dm3) || return zero(SqrtRational{BigInt})
    A = one(PFRational{Int16})
    sn = one(BigInt)
    _pf_3j_impl!(A, sn, dj1, dj2, dj3, dm1, dm2, dm3)
    iszero(sn) && return zero(SqrtRational{BigInt})
    sd = one(BigInt)
    rn = one(BigInt)
    rd = one(BigInt)
    extract_sr!(sn, sd, rn, rd, A)
    return SqrtRational(Base.unsafe_rational(sn, sd), Base.unsafe_rational(rn, rd))
end

function _pf_f3j(dj1::Int, dj2::Int, dj3::Int, dm1::Int, dm2::Int, dm3::Int)::Float64
    check_3j(dj1, dj2, dj3, dm1, dm2, dm3) || return zero(Float64)
    A = one(PFRational{Int16})
    sn = one(BigInt)
    _pf_3j_impl!(A, sn, dj1, dj2, dj3, dm1, dm2, dm3)
    iszero(sn) && return zero(Float64)
    sd = one(BigInt)
    rn = one(BigInt)
    rd = one(BigInt)
    extract_sr!(sn, sd, rn, rd, A)
    ans = BigFloat(sn) / BigFloat(sd) * sqrt(BigFloat(rn) / BigFloat(rd))
    return convert(Float64, ans)
end

function _pf_CG0(j1::Int, j2::Int, j3::Int)
    check_couple(2j1, 2j2, 2j3) || return zero(SqrtRational{BigInt})
    J::Int = j1 + j2 + j3
    isodd(J) && return zero(SqrtRational{BigInt})
    g = div(J, 2)
    Jsize::Int = max_prime_index(J + 1)
    A = pf_alloc(Int16, Jsize)
    _copy!(A, _pf_bin(g, j3))
    _mul!(A, _pf_bin(j3, g - j1))
    square!(A)
    _div!(A, _pf_bin(J + 1, 2j3 + 1))
    _div!(A, _pf_bin(2j3, J - 2j1))
    sn = one(BigInt)
    if isodd(g - j3)
        Base.GMP.MPZ.neg!(sn)
    end
    sd = one(BigInt)
    rn = one(BigInt)
    rd = one(BigInt)
    extract_sr!(sn, sd, rn, rd, A)
    return SqrtRational(Base.unsafe_rational(sn, sd), Base.unsafe_rational(rn, rd))
end

function _pf_fCG0(j1::Int, j2::Int, j3::Int)::Float64
    check_couple(2j1, 2j2, 2j3) || return zero(Float64)
    J::Int = j1 + j2 + j3
    isodd(J) && return zero(Float64)
    g = div(J, 2)
    Jsize::Int = max_prime_index(J + 1)
    A = pf_alloc(Int16, Jsize)
    _copy!(A, _pf_bin(g, j3))
    _mul!(A, _pf_bin(j3, g - j1))
    square!(A)
    _div!(A, _pf_bin(J + 1, 2j3 + 1))
    _div!(A, _pf_bin(2j3, J - 2j1))
    sn = one(BigInt)
    if isodd(g - j3)
        Base.GMP.MPZ.neg!(sn)
    end
    sd = one(BigInt)
    rn = one(BigInt)
    rd = one(BigInt)
    extract_sr!(sn, sd, rn, rd, A)
    ans = BigFloat(sn) / BigFloat(sd) * sqrt(BigFloat(rn) / BigFloat(rd))
    return convert(Float64, ans)
end

function _pf_6j(dj1::Int, dj2::Int, dj3::Int, dj4::Int, dj5::Int, dj6::Int)
    check_6j(dj1, dj2, dj3, dj4, dj5, dj6) || return zero(SqrtRational{BigInt})
    A = one(PFRational{Int16})
    sn = one(BigInt)
    _pf_6j_impl!(A, sn, dj1, dj2, dj3, dj4, dj5, dj6)
    iszero(sn) && return zero(SqrtRational{BigInt})
    sd = one(BigInt)
    rn = one(BigInt)
    rd = one(BigInt)
    extract_sr!(sn, sd, rn, rd, A)
    return SqrtRational(Base.unsafe_rational(sn, sd), Base.unsafe_rational(rn, rd))
end

function _pf_f6j(dj1::Int, dj2::Int, dj3::Int, dj4::Int, dj5::Int, dj6::Int)::Float64
    check_6j(dj1, dj2, dj3, dj4, dj5, dj6) || return zero(Float64)
    A = one(PFRational{Int16})
    sn = one(BigInt)
    _pf_6j_impl!(A, sn, dj1, dj2, dj3, dj4, dj5, dj6)
    iszero(sn) && return zero(Float64)
    sd = one(BigInt)
    rn = one(BigInt)
    rd = one(BigInt)
    extract_sr!(sn, sd, rn, rd, A)
    ans = BigFloat(sn) / BigFloat(sd) * sqrt(BigFloat(rn) / BigFloat(rd))
    return convert(Float64, ans)
end

"""
    eCG(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer)::SqrtRational{BigInt}
Exact CG coefficient using prime factorization algorithm.
"""
@inline function eCG(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer)::SqrtRational{BigInt}
    _pf_CG(convert(Int, dj1), convert(Int, dj2), convert(Int, dj3), convert(Int, dm1), convert(Int, dm2), convert(Int, dm3))
end

"""
    efCG(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer)::Float64
Exact CG coefficient (Float64) using prime factorization algorithm.
"""
@inline function efCG(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer)::Float64
    _pf_fCG(convert(Int, dj1), convert(Int, dj2), convert(Int, dj3), convert(Int, dm1), convert(Int, dm2), convert(Int, dm3))
end

"""
    e3j(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer)::SqrtRational{BigInt}
Exact 3j symbol using prime factorization algorithm.
"""
@inline function e3j(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer)::SqrtRational{BigInt}
    _pf_3j(convert(Int, dj1), convert(Int, dj2), convert(Int, dj3), convert(Int, dm1), convert(Int, dm2), convert(Int, dm3))
end

"""
    ef3j(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer)::Float64
Exact 3j symbol (Float64) using prime factorization algorithm.
"""
@inline function ef3j(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer)::Float64
    _pf_f3j(convert(Int, dj1), convert(Int, dj2), convert(Int, dj3), convert(Int, dm1), convert(Int, dm2), convert(Int, dm3))
end

"""
    eCG0(j1::Integer, j2::Integer, j3::Integer)::SqrtRational{BigInt}
Exact CG coefficient (`m1 = m2 = m3 = 0`) using prime factorization algorithm.
"""
@inline function eCG0(j1::Integer, j2::Integer, j3::Integer)::SqrtRational{BigInt}
    _pf_CG0(convert(Int, j1), convert(Int, j2), convert(Int, j3))
end

"""
    efCG0(j1::Integer, j2::Integer, j3::Integer)::Float64
Exact CG coefficient (`m1 = m2 = m3 = 0`) (Float64) using prime factorization algorithm.
"""
@inline function efCG0(j1::Integer, j2::Integer, j3::Integer)::Float64
    _pf_fCG0(convert(Int, j1), convert(Int, j2), convert(Int, j3))
end

"""
    e6j(j1::Integer, j2::Integer, j3::Integer, j4::Integer, j5::Integer, j6::Integer)::SqrtRational{BigInt}
Exact 6j symbol using prime factorization algorithm.
"""
@inline function e6j(j1::Integer, j2::Integer, j3::Integer, j4::Integer, j5::Integer, j6::Integer)::SqrtRational{BigInt}
    _pf_6j(convert(Int, j1), convert(Int, j2), convert(Int, j3), convert(Int, j4), convert(Int, j5), convert(Int, j6))
end

"""
    ef6j(j1::Integer, j2::Integer, j3::Integer, j4::Integer, j5::Integer, j6::Integer)::Float64
Exact 6j symbol (Float64) using prime factorization algorithm.
"""
@inline function ef6j(j1::Integer, j2::Integer, j3::Integer, j4::Integer, j5::Integer, j6::Integer)::Float64
    _pf_f6j(convert(Int, j1), convert(Int, j2), convert(Int, j3), convert(Int, j4), convert(Int, j5), convert(Int, j6))
end
