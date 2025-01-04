using Primes
using Base.GMP.MPZ

const c_tdiv_q_ui = MPZ.gmpz(:tdiv_q_ui)
const c_ui_pow_ui = MPZ.gmpz(:ui_pow_ui)
const c_divexact = MPZ.gmpz(:divexact)
const c_bin_uiui = MPZ.gmpz(:bin_uiui)
const c_addmul = MPZ.gmpz(:addmul)
@inline function __gmpz_tdiv_q_ui!(q::BigInt, n::BigInt, d::Culong)
    ccall(c_tdiv_q_ui, Culong, (Ref{BigInt}, Ref{BigInt}, Culong), q, n, d)
end
@inline function __gmpz_ui_pow_ui!(r::BigInt, b::Culong, e::Culong)
    ccall(c_ui_pow_ui, Cvoid, (Ref{BigInt}, Culong, Culong), r, b, e)
end
@inline function __gmpz_divexact!(q::BigInt, n::BigInt, d::BigInt)
    ccall(c_divexact, Cvoid, (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}), q, n, d)
end
@inline function __gmpz_bin_uiui!(r::BigInt, n::Culong, k::Culong)
    ccall(c_bin_uiui, Cvoid, (Ref{BigInt}, Culong, Culong), r, n, k)
end
@inline function __gmpz_addmul!(r::BigInt, a::BigInt, b::BigInt)
    ccall(c_addmul, Cvoid, (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}), r, a, b)
end

@inline function _bigbin(r::BigInt, n::Int, k::Int)
    __gmpz_bin_uiui!(r, Culong(n), Culong(k))
    return r
end

@inline _bigbin(n::Int, k::Int) = _bigbin(BigInt(), n, k)

@inline function _divgcd!(g::BigInt, n::BigInt, d::BigInt)
    MPZ.gcd!(g, n, d)
    __gmpz_divexact!(n, n, g)
    __gmpz_divexact!(d, d, g)
    return g
end

# t is buffer, modefiy `na, da` in place
@inline function __gmpq_add!(t::BigInt, na::BigInt, da::BigInt, nb::BigInt, db::BigInt)
    MPZ.mul!(na, db)
    __gmpz_addmul!(na, da, nb)
    MPZ.mul!(da, db)
    _divgcd!(t, na, da)
    return
end

@inline function __gmpq_add!(t::BigInt, a::Rational{BigInt}, nu::BigInt, de::BigInt)
    __gmpq_add!(t, a.num, a.den, nu, de)
end

@inline function __gmpq_add!(t::BigInt, a::Rational{BigInt}, b::Rational{BigInt})
    __gmpq_add!(t, a.num, a.den, b.num, b.den)
end

# simplify `x√t`, move the square factors of `t` to `x`
# assume `t` is positive
function simplify!(x::BigInt, t::BigInt, hint::Int)
    q = zero(BigInt)
    tt = one(BigInt)
    for ip in primes(hint)
        p = convert(Culong, ip)
        p > t && break
        r = __gmpz_tdiv_q_ui!(q, t, p)
        i::Culong = zero(Culong)
        while iszero(r)
            i += 1
            MPZ.set!(t, q)
            r = __gmpz_tdiv_q_ui!(q, t, p)
        end
        iszero(i) && continue
        ti, xi = divrem(i, Culong(2))
        xi == 1 && MPZ.mul_ui!(tt, p)
        if ti != 0
            __gmpz_ui_pow_ui!(q, p, ti)
            MPZ.mul!(x, q)
        end
    end
    MPZ.mul!(t, tt)
    return
end

# simplify s*√(n/d) -> s/g*√(n/d)
function simplify3!(g::BigInt, s::BigInt, n::BigInt, d::BigInt, hint::Int)
    _divgcd!(g, n, d)
    _divgcd!(g, s, d)
    MPZ.mul!(n, g)
    _divgcd!(g, n, d)
    MPZ.set_ui!(g, 1)
    simplify!(s, n, hint)
    simplify!(g, d, hint)
    return
end

# simplify sn/sd*√(rn/rd) -> sn/sd*√(rn/rd), t is buffer
function simplify4!(t::BigInt, sn::BigInt, sd::BigInt, rn::BigInt, rd::BigInt, hint::Int)
    _divgcd!(t, sn, sd)
    _divgcd!(t, rn, rd)
    simplify!(sn, rn, hint)
    simplify!(sd, rd, hint)
    _divgcd!(t, sn, sd)
    _divgcd!(t, sn, rd)
    MPZ.mul!(rn, t)
    _divgcd!(t, sd, rn)
    MPZ.mul!(rd, t)
    _divgcd!(t, rn, rd)
    return
end
