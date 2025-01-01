import Base: +, -, *, /, ==
using Markdown
using Primes

"""
    SqrtRational <: Real
Store exact value of a `Rational`'s square root.
It is stored in `s√r` format, and we do not simplify it during arithmetics.
You can use the `simplify` function to simplify it.
"""
struct SqrtRational{T} <: Real
    s::Rational{T}
    r::Rational{T}
    SqrtRational{T}(s::Rational{T}, r::Rational{T}) where {T} = begin
        r < 0 && throw(ArgumentError("invalid SqrtRational $s√$r"))
        (iszero(s) || iszero(r)) && return new(zero(s), one(r))
        new(s, r)
    end
end

SqrtRational(s::Rational{T1}, r::Rational{T2}) where {T1,T2} = SqrtRational{promote_type(T1, T2)}(promote(s, r)...)
SqrtRational(s::Integer, r::Rational{T}) where {T} = SqrtRational(promote(s, r)...)
SqrtRational(s::Rational{T}, r::Integer) where {T} = SqrtRational(promote(s, r)...)
SqrtRational(s::Integer, r::Integer) = SqrtRational(promote(Rational(s), r)...)
SqrtRational(s::Union{Integer,Rational}) = SqrtRational(s, one(s))

"""
    exact_sqrt(r::Union{Integer, Rational})
Get exact `√r` using `SqrtRational` type.
"""
exact_sqrt(r::Union{Integer,Rational}) = SqrtRational(one(r), r)

# some basic functions
Base.zero(::Type{SqrtRational{T}}) where {T} = SqrtRational(zero(T), one(T))
Base.one(::Type{SqrtRational{T}}) where {T} = SqrtRational(one(T), one(T))
Base.zero(x::SqrtRational) = zero(typeof(x))
Base.one(x::SqrtRational) = one(typeof(x))
Base.sign(x::SqrtRational) = sign(x.s)
Base.signbit(x::SqrtRational) = signbit(x.s)
Base.inv(x::SqrtRational) = SqrtRational(inv(x.s), inv(x.r))

# override some operators
+(x::SqrtRational) = x
-(x::SqrtRational) = SqrtRational(-x.s, x.r)
*(x::SqrtRational, y::SqrtRational) = SqrtRational(x.s * y.s, x.r * y.r)
*(x::SqrtRational, y::Union{Integer,Rational}) = SqrtRational(x.s * y, x.r)
*(x::Union{Integer,Rational}, y::SqrtRational) = y * x
/(x::SqrtRational, y::SqrtRational) = SqrtRational(x.s / y.s, x.r / y.r)
/(x::SqrtRational, y::Union{Integer,Rational}) = SqrtRational(x.s / y, x.r)
/(x::Union{Integer,Rational}, y::SqrtRational) = SqrtRational(x / y.s, inv(y.r))
==(x::SqrtRational, y::SqrtRational) = begin
    if y.s == 0
        return x.s == 0
    else
        isone((x.s / y.s)^2 * (x.r / y.r))
    end
end
==(x::SqrtRational, y::Union{Integer,Rational}) = (x == SqrtRational(y))
==(x::Union{Integer,Rational}, y::SqrtRational) = (SqrtRational(x) == y)

# Only if `r == 1` in `s√r`, the add operator can work.
# We do not simplify `s√r` every time, but in this function,
# we need to simplify it first, and then do the `+` operator.
+(x::SqrtRational, y::Union{Integer,Rational}) = begin
    r = x.r
    t1 = isqrt(numerator(r))
    t1 * t1 != numerator(r) && throw(ArgumentError("cannot simplify $(x)"))
    t2 = isqrt(denominator(r))
    t2 * t2 != denominator(r) && throw(ArgumentError("cannot simplify $(x)"))
    return x.s * (t1 // t2) + y
end
+(x::Union{Integer,Rational}, y::SqrtRational) = y + x
+(x::SqrtRational{T1}, y::SqrtRational{T2}) where {T1,T2} = begin
    if x == 0
        return y
    else
        return x * (one(promote_type(T1, T2)) + y / x)
    end
end
-(x::SqrtRational, y::Union{Integer,Rational}) = (x + (-y))
-(x::Union{Integer,Rational}, y::SqrtRational) = (x + (-y))
-(x::SqrtRational, y::SqrtRational) = (x + (-y))

# widen
Base.widen(x::SqrtRational) = SqrtRational(widen(x.s), widen(x.r))

# convert
"""
    float(x::SqrtRational)::BigFloat
Convert a `SqrtRational` to `BigFloat`.
"""
Base.float(x::SqrtRational) = convert(BigFloat, x.s) * sqrt(convert(BigFloat, x.r))

# show
Base.show(io::IO, ::MIME"text/plain", x::SqrtRational) = begin
    if isone(x.r)
        if denominator(x.s) == 1
            print(io, numerator(x.s))
        else
            print(io, x.s)
        end
        return
    elseif isone(x.s)
        if denominator(x.r) == 1
            print(io, "√$(numerator(x.r))")
        else
            print(io, "√($(x.r))")
        end
        return
    elseif isone(-x.s)
        if denominator(x.r) == 1
            print(io, "-√$(numerator(x.r))")
        else
            print(io, "-√($(x.r))")
        end
        return
    end
    to_show::String = ""
    if denominator(x.s) == 1
        to_show = string(numerator(x.s))
    else
        to_show = "$(x.s)"
    end
    if denominator(x.r) == 1
        to_show *= "√$(numerator(x.r))"
    else
        to_show *= "√($(x.r))"
    end
    print(io, to_show)
end

Base.show(io::IO, x::SqrtRational) = show(io::IO, "text/plain", x)

Base.show(io::IO, ::MIME"text/markdown", x::SqrtRational) = begin
    if isone(x.r)
        if denominator(x.s) == 1
            show(io, "text/markdown", Markdown.parse("``$(numerator(x.s))``"))
        else
            show(io, "text/markdown", Markdown.parse("``\\frac{$(numerator(x.s))}{$(denominator(x.s))}``"))
        end
        return
    elseif isone(x.s)
        if denominator(x.r) == 1
            show(io, "text/markdown", Markdown.parse("``\\sqrt{$(numerator(x.r))}``"))
        else
            show(io, "text/markdown", Markdown.parse("``\\sqrt{\\frac{$(numerator(x.r))}{$(denominator(x.r))}}``"))
        end
        return
    elseif isone(-x.s)
        if denominator(x.r) == 1
            show(io, "text/markdown", Markdown.parse("``-\\sqrt{$(numerator(x.r))}``"))
        else
            show(io, "text/markdown", Markdown.parse("``-\\sqrt{\\frac{$(numerator(x.r))}{$(denominator(x.r))}}``"))
        end
        return
    end
    to_show::String = ""
    if denominator(x.s) == 1
        to_show = string(numerator(x.s))
    else
        to_show = "\\frac{$(numerator(x.s))}{$(denominator(x.s))}"
    end
    if denominator(x.r) == 1
        to_show *= "\\sqrt{$(numerator(x.r))}"
    else
        to_show *= "\\sqrt{\\frac{$(numerator(x.r))}{$(denominator(x.r))}}"
    end
    show(io, "text/markdown", Markdown.parse("``$to_show``"))
end

const c_tdiv_q_ui = Base.GMP.MPZ.gmpz(:tdiv_q_ui)
@inline function __gmpz_tdiv_q_ui!(q::BigInt, n::BigInt, d::Culong)
    ccall(c_tdiv_q_ui, Culong, (Ref{BigInt}, Ref{BigInt}, Culong), q, n, d)
end

function simplify!(n::BigInt, hint::Int)
    iszero(n) && return zero(BigInt), zero(BigInt)
    s = sign(n)
    s < 0 && Base.GMP.MPZ.neg!(n)
    isone(n) && return s, n
    if hint == 0
        x = one(BigInt)
        t = one(BigInt)
        for (f, i) in factor(n)
            ti, xi = divrem(i, 2)
            xi == 1 && (x = x * f)
            ti != 0 && (t = t * f^ti)
        end
        return s * x, t
    end
    x = one(BigInt)
    t = one(BigInt)
    q = zero(BigInt)
    for p in primes(hint)
        p = convert(Culong, p)
        if p > n
            break
        end
        r = __gmpz_tdiv_q_ui!(q, n, p)
        i = zero(Culong)
        while r == 0
            i += 1
            Base.GMP.MPZ.set!(n, q)
            r = __gmpz_tdiv_q_ui!(q, n, p)
        end
        ti, xi = divrem(i, 2)
        xi == 1 && Base.GMP.MPZ.mul_ui!(x, p)
        ti != 0 && Base.GMP.MPZ.mul!(t, big(p)^ti)
    end
    Base.GMP.MPZ.mul!(x, n)
    s < 0 && Base.GMP.MPZ.neg!(x)
    return s * x, t
end

function simplify!(x::SqrtRational{BigInt}, hint::Int)
    nx, nt = simplify!(numerator(x.r), hint)
    dx, dt = simplify!(denominator(x.r), hint)
    nt *= numerator(x.s)
    dt *= denominator(x.s)
    a0 = gcd(nt, dt)
    nt = div(nt, a0)
    dt = div(dt, a0)
    a1 = gcd(nx, dt)
    a2 = gcd(dx, nt)
    dt = div(dt, a1)
    nt = div(nt, a2)
    dx = div(dx, a2) * a1
    nx = div(nx, a1) * a2
    return SqrtRational(nt // dt, nx // dx)
end

"""
    simplify(n::Integer)
Simplify a integer `n = x * t^2` to `(x, t)`
"""
function simplify(n::Integer)
    s = sign(n)
    n = s * n
    x = one(n)
    t = one(n)
    for (f, i) in factor(n)
        ti, xi = divrem(i, 2)
        xi == 1 && (x = x * f)
        ti != 0 && (t = t * f^ti)
    end
    return s * x, t
end

"""
    simplify(x::SqrtRational)
Simplify a SqrtRational.
"""
function simplify(x::SqrtRational)
    nx, nt = simplify(numerator(x.r))
    dx, dt = simplify(denominator(x.r))
    nt *= numerator(x.s)
    dt *= denominator(x.s)
    a0 = gcd(nt, dt)
    nt = div(nt, a0)
    dt = div(dt, a0)
    a1 = gcd(nx, dt)
    a2 = gcd(dx, nt)
    dt = div(dt, a1)
    nt = div(nt, a2)
    dx = div(dx, a2) * a1
    nx = div(nx, a1) * a2
    return SqrtRational(nt // dt, nx // dx)
end