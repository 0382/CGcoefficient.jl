import Base:+,-,*,/,==
using Markdown

"""
    SqrtRational <: Real
A type to store exact value of square root of a `Rational`: ±√p. To store infinity values, use `Rational{BigInt}` to represent the `p`.
Basic functions such as `zero, one, inv, *, /` are supported for this type.
"""
struct SqrtRational <: Real
    sign::Bool
    value::Rational{BigInt}
    SqrtRational(s::Bool, x::Rational{BigInt}) = begin
        isinf(x) || x < 0 && throw(ArgumentError("invalid SqrtRational  √($x)"))
        new(s, x)
    end
end

SqrtRational(s::Bool, x::Integer) = SqrtRational(s, Rational{BigInt}(big(x)))
SqrtRational(s::Bool, x::Rational) = SqrtRational(s, big(x))
SqrtRational(x::Integer) = SqrtRational(false, Rational{BigInt}(big(x)))
SqrtRational(x::Rational) = SqrtRational(false, Rational{BigInt}(big(x)))

# some basic functions
Base.zero(::Type{SqrtRational}) = SqrtRational(zero(Rational{BigInt}))
Base.zero(x::SqrtRational) = SqrtRational(zero(Rational{BigInt}))
Base.one(::Type{SqrtRational}) = SqrtRational(one(Rational{BigInt}))
Base.one(x::SqrtRational) = SqrtRational(one(Rational{BigInt}))
Base.sign(x::SqrtRational) = iszero(x.value) ? 0 : x.sign ? -1 : 1
Base.signbit(x::SqrtRational) = x.sign
Base.inv(x::SqrtRational) = SqrtRational(x.sign, inv(x.value))

# override some operators
+(x::SqrtRational) = x
-(x::SqrtRational) = SqrtRational(!(x.sign), x.value)
*(x::SqrtRational, y::SqrtRational) = SqrtRational(x.sign ⊻ y.sign, x.value * y.value)
*(x::SqrtRational, y::Union{Rational, Integer}) = x * SqrtRational(signbit(y), abs(big(y))^2)
*(x::Union{Rational, Integer}, y::SqrtRational) = y * x
/(x::SqrtRational, y::Union{SqrtRational, Rational, Integer}) = x * inv(y)
/(x::Union{Rational, Integer}, y::SqrtRational) = SqrtRational(signbit(x), abs(big(x))^2) / y
==(x::SqrtRational, y::SqrtRational) = (sign(x) == sign(y)) & (x.value == y.value)

Base.float(x::SqrtRational) = sign(x) * sqrt(float(x.value))

"""
simplify a integer `n = x * t^2` to `(x, t)`
"""
function simplify(x::T) where T<:Integer
    s = sign(x)
    x = s * x
    out = one(T)
    i = 2
    while i <= floor(T, sqrt(x))
        i2 = i * i
        tp = div(x, i2)
        if tp * i2 == x
            out *= i
            x = tp
            i = 2
        else
            i += 1
        end
    end
    return out, s * x
end

"""
simplify a SqrtRational `√(p/q)` to x/t√(p1/q1)
"""
function simplify(x::SqrtRational)
    s = signbit(x)
    nu_out, nu_in = simplify(numerator(x.value))
    de_out, de_in = simplify(denominator(x.value))
    return s, nu_out, nu_in, de_out, de_in
end

Base.show(io::IO, ::MIME"text/plain", x::SqrtRational) = begin
    s, nu_out, nu_in, de_out, de_in = simplify(x)
    if nu_in == 0
        print(io, 0)
        return
    end
    s && print(io, '-')
    if nu_in == 1
        print(io, nu_out)
    elseif nu_out == 1
        print(io, "√", nu_in)
    else
        print(io, nu_out, "√", nu_in)
    end
    if de_out == de_in == 1
        return
    elseif de_out == 1
        print(io, "/√", de_in)
    elseif de_in == 1
        print(io, "/", de_out)
    else
        print(io, "/", de_out, "√", de_in)
    end
    return
end

Base.show(io::IO, x::SqrtRational) = show(io::IO, "text/plain", x)

Base.show(io::IO, ::MIME"text/markdown", x::SqrtRational) = begin
    s, nu_out, nu_in, de_out, de_in = simplify(x)
    if nu_in == 0
        show(io, "text/markdown", Markdown.parse("``0``"))
        return
    end
    to_show::String = ""
    if nu_in == 1
        to_show = "$nu_out"
    elseif nu_out == 1
        to_show = "\\sqrt{$nu_in}"
    else
        to_show = "$nu_out\\sqrt{$nu_in}"
    end
    if de_out == de_in == 1
        nothing
    elseif de_out == 1
        to_show = "\\dfrac{$to_show}{\\sqrt{$de_in}}"
    elseif de_in == 1
        to_show = "\\dfrac{$to_show}{$de_out}"
    else
        to_show = "\\dfrac{$to_show}{$de_out\\sqrt{$de_in}}"
    end
    to_show = (s ? "-" : "") * to_show
    show(io, "text/markdown", Markdown.parse("``$to_show``"))
end
