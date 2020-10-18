import Base:+,-,*,/,==
using Markdown

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

# 一些基础的函数
Base.zero(::Type{SqrtRational}) = SqrtRational(zero(Rational{BigInt}))
Base.zero(x::SqrtRational) = SqrtRational(zero(Rational{BigInt}))
Base.one(::Type{SqrtRational}) = SqrtRational(one(Rational{BigInt}))
Base.one(x::SqrtRational) = SqrtRational(one(Rational{BigInt}))
Base.sign(x::SqrtRational) = iszero(x.value) ? 0 : x.sign ? -1 : 1
Base.signbit(x::SqrtRational) = x.sign
Base.inv(x::SqrtRational) = SqrtRational(x.sign, inv(x.value))

# 运算符重载
+(x::SqrtRational) = x
-(x::SqrtRational) = SqrtRational(!(x.sign), x.value)
*(x::SqrtRational, y::SqrtRational) = SqrtRational(x.sign ⊻ y.sign, x.value * y.value)
*(x::SqrtRational, y::Union{Rational, Integer}) = x * SqrtRational(signbit(y), abs(big(y))^2)
*(x::Union{Rational, Integer}, y::SqrtRational) = y * x
/(x::SqrtRational, y::Union{SqrtRational, Rational, Integer}) = x * inv(y)
/(x::Union{Rational, Integer}, y::SqrtRational) = SqrtRational(signbit(x), abs(big(x))^2) / y
==(x::SqrtRational, y::SqrtRational) = (sign(x) == sign(y)) & (x.value == y.value)

Base.float(x::SqrtRational) = sign(x) * sqrt(float(x.value))

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

function simplify(x::SqrtRational)
    s = signbit(x)
    nu_o, nu_i = simplify(numerator(x.value))
    de_o, de_i = simplify(denominator(x.value))
    return s, nu_o, nu_i, de_o, de_i
end

Base.show(io::IO, ::MIME"text/plain", x::SqrtRational) = begin
    s, nu_o, nu_i, de_o, de_i = simplify(x)
    if nu_i == 0
        print(io, 0)
        return
    end
    s && print(io, '-')
    if nu_i == 1
        print(io, nu_o)
    elseif nu_o == 1
        print(io, "√", nu_i)
    else
        print(io, nu_o, "√", nu_i)
    end
    if de_o == de_i == 1
        return
    elseif de_o == 1
        print(io, "/√", de_i)
    elseif de_i == 1
        print(io, "/", de_o)
    else
        print(io, "/", de_o, "√", de_i)
    end
    return
end

Base.show(io::IO, x::SqrtRational) = show(io::IO, "text/plain", x)

Base.show(io::IO, ::MIME"text/markdown", x::SqrtRational) = begin
    s, nu_o, nu_i, de_o, de_i = simplify(x)
    if nu_i == 0
        show(io, "text/markdown", Markdown.parse("``0``"))
        return
    end
    _to_show::String = ""
    if nu_i == 1
        _to_show = "$nu_o"
    elseif nu_o == 1
        _to_show = "\\sqrt{$nu_i}"
    else
        _to_show = "$nu_o\\sqrt{$nu_i}"
    end
    if de_o == de_i == 1
        nothing
    elseif de_o == 1
        _to_show = "\\dfrac{$_to_show}{\\sqrt{$de_i}}"
    elseif de_i == 1
        _to_show = "\\dfrac{$_to_show}{$de_o}"
    else
        _to_show = "\\dfrac{$_to_show}{$de_o\\sqrt{$de_i}}"
    end
    _to_show = (s ? "-" : "") * _to_show
    show(io, "text/markdown", Markdown.parse("``$_to_show``"))
end
