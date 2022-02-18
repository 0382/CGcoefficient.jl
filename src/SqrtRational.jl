import Base:+,-,*,/,==
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
    SqrtRational{T}(s::Rational{T}, r::Rational{T}) where T = begin
        r <= 0 && throw(ArgumentError("invalid SqrtRational $s√$r"))
        new(s, r)
    end
end

SqrtRational(s::Rational{T1}, r::Rational{T2}) where {T1, T2} = SqrtRational{promote_type(T1, T2)}(promote(s, r)...)
SqrtRational(s::Integer, r::Rational{T}) where T = SqrtRational(promote(s, r)...)
SqrtRational(s::Rational{T}, r::Integer) where T = SqrtRational(promote(s, r)...)
SqrtRational(s::Integer, r::Integer) = SqrtRational(promote(Rational(s), r)...)
SqrtRational(s::Union{Integer, Rational}) = SqrtRational(s, one(s))

"""
    exact_sqrt(r::Union{Integer, Rational})
Get exact `√r` using `SqrtRational` type.
"""
exact_sqrt(r::Union{Integer, Rational}) = SqrtRational(one(r), r)

# some basic functions
Base.zero(::Type{SqrtRational{T}}) where T = SqrtRational(zero(T), one(T))
Base.one(::Type{SqrtRational{T}}) where T = SqrtRational(one(T), one(T))
Base.zero(x::SqrtRational) = zero(typeof(x))
Base.one(x::SqrtRational) = one(typeof(x))
Base.sign(x::SqrtRational) = sign(x.s)
Base.signbit(x::SqrtRational) = signbit(x.s)
Base.inv(x::SqrtRational) = SqrtRational(inv(x.s), inv(x.r))

# override some operators
+(x::SqrtRational) = x
-(x::SqrtRational) = SqrtRational(-x.s, x.r)
*(x::SqrtRational, y::SqrtRational) = SqrtRational(x.s * y.s, x.r * y.r)
*(x::SqrtRational, y::Union{Integer, Rational}) = SqrtRational(x.s * y, x.r)
*(x::Union{Integer, Rational}, y::SqrtRational) = y * x
/(x::SqrtRational, y::SqrtRational) = SqrtRational(x.s/y.s, x.r/y.r)
/(x::SqrtRational, y::Union{Integer, Rational}) = SqrtRational(x.s/y, x.r)
/(x::Union{Integer, Rational}, y::SqrtRational) = SqrtRational(x/y.s, inv(y.r))
==(x::SqrtRational, y::SqrtRational) = begin
    if y.s == 0
        return x.s == 0
    else
        isone((x.s/y.s)^2 * (x.r/y.r))
    end
end
==(x::SqrtRational, y::Union{Integer, Rational}) = (x == SqrtRational(y))
==(x::Union{Integer, Rational}, y::SqrtRational) = (SqrtRational(x) == y)

"""
Only if `r == 1` in `s√r`, the add operator can work.
We do not simplify `s√r` every time, but in this function,
we need to simplify it first, and then do the `+` operator.
"""
+(x::SqrtRational, y::Union{Integer, Rational}) = begin
    r = x.r
    t1 = isqrt(numerator(r))
    t1 * t1 != numerator(r) && throw(ArgumentError("cannot simplify $(x.s)√$(x.r)"))
    t2 = isqrt(denominator(r))
    t2 * t2 != denominator(r) && throw(ArgumentError("cannot simplify $(x.s)√$(x.r)"))
    return x.s * (t1//t2) + y
end
+(x::Union{Integer, Rational}, y::SqrtRational) = y + x
+(x::SqrtRational{T1}, y::SqrtRational{T2}) where {T1, T2} = begin
    if x == 0
        return y
    else
        return x * (one(promote_type(T1, T2)) + y/x)
    end
end
-(x::SqrtRational, y::Union{Integer, Rational}) = (x + (-y))
-(x::Union{Integer, Rational}, y::SqrtRational) = (x + (-y))
-(x::SqrtRational, y::SqrtRational) = (x + (-y))

# widen
Base.widen(x::SqrtRational) = SqrtRational(widen(x.s), widen(x.r))

# convert
"""
Although we use `BigInt` in calculation, we do not convert it to `BigFloat`.
Because this package is designed for numeric calculation, `BigFloat` is unnecessary.
"""
Base.float(x::SqrtRational) = Float64(x.s) * sqrt(Float64(x.r))

# show
Base.show(io::IO, ::MIME"text/plain", x::SqrtRational) = begin
    if iszero(x.s)
        print(io, 0)
        return
    elseif isone(x.s)
        if denominator(x.r) == 1
            print(io, "√$(numerator(x.r))")
        else
            print(io, "√($(x.r))")
        end
        return
    end
    to_show::String = ""
    if denominator(x.s) == 1
        to_show = string(numerator(x.s))
    else
        to_show = "($(x.s))"
    end
    if denominator(x.r) == 1
        if numerator(x.r) != 1
            to_show = "$to_show√$(numerator(x.r))"
        end
    else
        to_show = "$to_show√($(x.r))"
    end
    print(io, to_show)
end

Base.show(io::IO, x::SqrtRational) = show(io::IO, "text/plain", x)

Base.show(io::IO, ::MIME"text/markdown", x::SqrtRational) = begin
    if iszero(x.s)
        show(io, "text/markdown", Markdown.parse("``0``"))
        return
    elseif isone(x.s)
        if denominator(x.r) == 1
            show(io, "text/markdown", Markdown.parse("``\\sqrt{$(numerator(x.r))}``"))
        else
            show(io, "text/markdown", Markdown.parse("``\\sqrt{\\frac{$(numerator(x.r))}{$(denominator(x.r))}}``"))
        end
    end
    to_show::String = ""
    if denominator(x.s) == 1
        to_show = string(numerator(x.s))
    else
        to_show = "\\frac{$(numerator(x.s))}{$(denominator(x.s))}"
    end
    if denominator(x.r) == 1
        if numerator(x.r) != 1
            to_show = "$to_show\\sqrt{$(numerator(x.r))}"
        end
    else
        to_show = "$to_show\\sqrt{\\frac{$(numerator(x.r))}{$(denominator(x.r))}}"
    end
    show(io, "text/markdown", Markdown.parse("``$to_show``"))
end

"""
    simplify(n::Integer)
Simplify a integer `n = x * t^2` to `(x, t)`
"""
function simplify(n::T) where T<:Integer
    s = sign(n)
    n = s * n
    t = one(T)
    i = convert(T, 2)
    while i <= isqrt(n)
        i2 = i * i
        temp = div(n, i2)
        if temp * i2 == n
            t *= i
            n = temp
        else
            i = nextprime(i, 2)
        end
    end
    return s * n, t
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
    return SqrtRational(nt//dt, nx//dx)
end