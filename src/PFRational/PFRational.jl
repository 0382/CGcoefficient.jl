import Base:*,/,^

@doc raw"""
    PFRational{T<:Integer} <: Real
A positive rational number in prime factorization form.
```math
q = \prod_{i=1}^{n} p_i^{e_i}
```
where ``p_i`` is the ``i``-th prime number and ``e_i`` can be negative.
"""
struct PFRational{T<:Integer} <: Real
    e::Vector{T}
    PFRational{T}(e::Vector{T}) where {T} = new{T}(e)
end

PFRational(e::Vector{T}) where {T<:Integer} = PFRational{T}(copy(e))

Base.one(::Type{PFRational{T}}) where {T<:Integer} = PFRational{T}(zeros(T, 0))
Base.copy(x::PFRational) = PFRational(copy(x.e))
Base.eltype(::Type{PFRational{T}}) where {T<:Integer} = T

pf_alloc(T::Type{<:Integer}, n::Integer) = PFRational{T}(zeros(T, n))
pf_alloc!(x::PFRational{T}, n::Integer) where {T<:Integer} = begin
    resize!(x.e, n)
    fill!(x.e, T(0))
end

Base.inv(x::PFRational) = PFRational(-x.e)

*(x::PFRational, y::PFRational) = begin
    T = promote_type(eltype(x.e), eltype(y.e))
    e = zeros(T, max(length(x.e), length(y.e)))
    copyto!(e, x.e)
    for i in eachindex(y.e)
        e[i] += y.e[i]
    end
    return PFRational{T}(e)
end

/(x::PFRational, y::PFRational) = begin
    T = promote_type(eltype(x.e), eltype(y.e))
    e = zeros(T, max(length(x.e), length(y.e)))
    copyto!(e, x.e)
    for i in eachindex(y.e)
        e[i] -= y.e[i]
    end
    return PFRational{T}(e)
end

^(x::PFRational, n::Integer) = PFRational(x.e .* n)
square!(x::PFRational) = begin
    for i in eachindex(x.e)
        x.e[i] *= 2
    end
end

"""
    gcd(x::PFRational, y::PFRational)
gcd of rational means gcd of numerator and lcm of denominator
"""
Base.gcd(x::PFRational, y::PFRational) = begin
    T = promote_type(eltype(x.e), eltype(y.e))
    e = zeros(T, max(length(x.e), length(y.e)))
    copyto!(e, x.e)
    for i in eachindex(y.e)
        e[i] = min(e[i], y.e[i])
    end
    return PFRational{T}(e)
end

"""
    lcm(x::PFRational, y::PFRational)
lcm of rational means lcm of numerator and gcd of denominator
"""
Base.lcm(x::PFRational, y::PFRational) = begin
    T = promote_type(eltype(x.e), eltype(y.e))
    e = zeros(T, max(length(x.e), length(y.e)))
    copyto!(e, x.e)
    for i in eachindex(y.e)
        e[i] = max(e[i], y.e[i])
    end
    return PFRational{T}(e)
end

# assume length(y) <= length(x)
_copy!(x::PFRational{T}, y::PFRational{T}) where {T<:Integer} = begin
    for i in eachindex(y.e)
        x.e[i] = y.e[i]
    end
    for i in length(y.e)+1:length(x.e)
        x.e[i] = 0
    end
end

# assume length(y) <= length(x)
_mul!(x::PFRational{T}, y::PFRational{T}) where {T<:Integer} = begin
    for i in eachindex(y.e)
        x.e[i] += y.e[i]
    end
end

# assume length(y) <= length(x)
_div!(x::PFRational{T}, y::PFRational{T}) where {T<:Integer} = begin
    for i in eachindex(y.e)
        x.e[i] -= y.e[i]
    end
end

# assume length(y) <= length(x)
_gcd!(x::PFRational{T}, y::PFRational{T}) where {T<:Integer} = begin
    for i in eachindex(y.e)
        x.e[i] = min(x.e[i], y.e[i])
    end
    for i in length(y.e)+1:length(x.e)
        x.e[i] = min(x.e[i], 0)
    end
end

# assume length(y) <= length(x)
_lcm!(x::PFRational{T}, y::PFRational{T}) where {T<:Integer} = begin
    for i in eachindex(y.e)
        x.e[i] = max(x.e[i], y.e[i])
    end
    for i in length(y.e)+1:length(x.e)
        x.e[i] = max(x.e[i], 0)
    end
end

# positive integers in prime factorization form
_pf_integers = PFRational{Int16}[one(PFRational{Int16})]
# binomial coefficients in prime factorization form
_pf_bin_data = PFRational{Int16}[one(PFRational{Int16})]
# maximum `n` for binomial coefficients
_pf_bin_nmax = Int(0)

function _extend_pf_integers(n::Integer)
    n <= length(_pf_integers) && return
    _extend_prime_table_to(n)
    sizehint!(_pf_integers, n)
    maxsize = length(_prime_table)
    e = zeros(Int16, maxsize)
    for i in length(_pf_integers)+1:n
        fill!(e, Int16(0))
        idx = 0
        t = i
        while t > 1
            idx += 1
            p = nth_stored_prime(idx)
            while t % p == 0
                e[idx] += 1
                t = div(t, p)
            end
        end
        push!(_pf_integers, PFRational{Int16}(deepcopy(e[1:idx])))
    end
end

@inline pf_int(n::Integer)::PFRational{Int16} = _pf_integers[n]::PFRational{Int16}
@inline pf_bin_nmax()::Int = _pf_bin_nmax::Int
@inline pf_bin_data()::Vector{PFRational{Int16}} = _pf_bin_data::Vector{PFRational{Int16}}
@inline function _pf_bin(n::Integer, k::Integer)::PFRational{Int16}
    k = min(k, n - k)
    return @inbounds pf_bin_data()[binomial_index(n, k)]
end

"""
    pf_binomial(n::Integer, k::Integer)::PFRational{Int16}
Return the stored binomial coefficient `C(n, k)` in prime factorization form.
"""
function pf_binomial(n::Integer, k::Integer)::PFRational{Int16}
    0 <= k <= n <= pf_bin_nmax() || throw(ArgumentError("invalid binomial arguments"))
    k = min(k, n - k)
    return pf_bin_data()[binomial_index(n, k)]
end


function _extend_pf_binomial(n::Integer)
    n <= pf_bin_nmax() && return
    _extend_prime_table_to(n)
    _extend_pf_integers(n)
    sizehint!(_pf_bin_data, binomial_data_size(n))
    pf = pf_alloc(Int16, length(_prime_table))
    for i in pf_bin_nmax()+1:n
        fill!(pf.e, Int16(0))
        push!(_pf_bin_data, one(PFRational{Int16}))
        for j in 1:div(i, 2)
            _mul!(pf, pf_int(i - j + 1))
            _div!(pf, pf_int(j))
            size = length(pf.e)
            while size > 0 && pf.e[size] == 0
                size -= 1
            end
            push!(_pf_bin_data, PFRational{Int16}(deepcopy(pf.e[1:size])))
        end
    end
    global _pf_bin_nmax = n
end

_numerator!(num::BigInt, x::PFRational) = begin
    MPZ.set_ui!(num, 1)
    for i in eachindex(x.e)
        if x.e[i] > 0
            mul_p_pow_k!(num, i, convert(Int, x.e[i]))
        end
    end
    return num
end

_denominator!(den::BigInt, x::PFRational) = begin
    MPZ.set_ui!(den, 1)
    for i in eachindex(x.e)
        if x.e[i] < 0
            mul_p_pow_k!(den, i, -convert(Int, x.e[i]))
        end
    end
    return den
end

_rational!(num::BigInt, den::BigInt, x::PFRational) = begin
    MPZ.set_ui!(num, 1)
    MPZ.set_ui!(den, 1)
    for i in eachindex(x.e)
        if x.e[i] > 0
            mul_p_pow_k!(num, i, convert(Int, x.e[i]))
        elseif x.e[i] < 0
            mul_p_pow_k!(den, i, -convert(Int, x.e[i]))
        end
    end
end

Base.numerator(x::PFRational) = _numerator!(BigInt(0), x)
Base.denominator(x::PFRational) = _denominator!(BigInt(0), x)

Base.convert(::Type{Rational}, x::PFRational) = begin
    num = BigInt(0)
    den = BigInt(0)
    _rational!(num, den, x)
    return Base.unsafe_rational(num, den)
end

Base.show(io::IO, x::PFRational{T}) where T = begin
    _extend_prime_table(length(x.e))
    print(io, "PFRational{", T, "}(")
    first = true
    for i in eachindex(x.e)
        if x.e[i] != 0
            p = nth_stored_prime(i)
            if !first
                print(io, " * ")
            end
            first = false
            print(io, p)
            if x.e[i] != 1
                print(io, "^", x.e[i])
            end
        end
    end
    print(io, ")")
end

"""
    wigner_init_pf(n::Integer, mode::AbstractString, rank::Integer)
Initialize the prime factorization tables for Wigner symbols.
The parameters are similar to `wigner_inif_float`.
You must call this function before using the `e` and `ef` versions of Wigner symbols.
"""
function wigner_init_pf(n::Integer, mode::AbstractString, rank::Integer)
    nmax = 0
    if mode == "Jmax"
        rank == 3 && (nmax = 3n + 1)
        rank == 6 && (nmax = 4n + 1)
        rank == 9 && (nmax = 5n + 1)
        (rank != 3 && rank != 6 && rank != 9) && throw(ArgumentError("invalid rank: $rank"))
    elseif mode == "2bjmax"
        rank == 3 && (nmax = 2n + 1)
        rank == 6 && (nmax = 3n + 1)
        rank == 9 && (nmax = 4n + 1)
        (rank != 3 && rank != 6 && rank != 9) && throw(ArgumentError("invalid rank: $rank"))
    elseif mode == "nmax"
        nmax = n
        rank = 3
    else
        throw(ArgumentError("invalid mode: $mode"))
    end
    # this is just a rough guess, if it throws an error, we can increase it
    approx_times = Dict(3=>8, 6=>10, 9=>12)
    _extend_prime_table_to(nmax)
    size = length(_prime_table)
    _extend_prime_pow_table(size, nmax, get(approx_times, rank, 8))
    _extend_pf_integers(nmax)
    _extend_pf_binomial(nmax)
    nothing
end
