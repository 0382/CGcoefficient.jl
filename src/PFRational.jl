# import Base:*,/,^
# using StaticArrays

# struct PFRational{T<:Integer} <: Real
#     e::Vector{T}
#     PFRational{T}(e::Vector{T}) where {T} = new{T}(e)
# end

# PFRational(e::Vector{T}) where {T<:Integer} = PFRational{T}(copy(e))

# Base.one(::Type{PFRational{T}}) where {T<:Integer} = PFRational{T}(Vector{T}(undef, 0))
# Base.copy(x::PFRational) = PFRational(copy(x.e))

# pf_alloc(T::Type{<:Integer}, n::Integer) = PFRational{T}(Vector{T}(undef, n))

# Base.inv(x::PFRational) = begin
#     e = copy(x.e)
#     for i in eachindex(e)
#         e[i] = -e[i]
#     end
#     return PFRational(e)
# end

# *(x::PFRational, y::PFRational) = begin
#     T = promote_type(eltype(x.e), eltype(y.e))
#     e = zeros(T, max(length(x.e), length(y.e)))
#     copyto!(e, x.e)
#     for i in eachindex(y.e)
#         e[i] += y.e[i]
#     end
#     return PFRational{T}(e)
# end

# /(x::PFRational, y::PFRational) = begin
#     T = promote_type(eltype(x.e), eltype(y.e))
#     e = zeros(T, max(length(x.e), length(y.e)))
#     copyto!(e, x.e)
#     for i in eachindex(y.e)
#         e[i] -= y.e[i]
#     end
#     return PFRational{T}(e)
# end

# ^(x::PFRational, n::Integer) = begin
#     e = copy(x.e)
#     for i in eachindex(e)
#         e[i] *= n
#     end
#     return PFRational(e)
# end

# """
#     gcd(x::PFRational, y::PFRational)
# gcd of rational means gcd of numerator and lcm of denominator
# """
# Base.gcd(x::PFRational, y::PFRational) = begin
#     T = promote_type(eltype(x.e), eltype(y.e))
#     e = zeros(T, max(length(x.e), length(y.e)))
#     copyto!(e, x.e)
#     for i in eachindex(y.e)
#         e[i] = min(e[i], y.e[i])
#     end
#     return PFRational{T}(e)
# end

# """
#     lcm(x::PFRational, y::PFRational)
# lcm of rational means lcm of numerator and gcd of denominator
# """
# Base.lcm(x::PFRational, y::PFRational) = begin
#     T = promote_type(eltype(x.e), eltype(y.e))
#     e = zeros(T, max(length(x.e), length(y.e)))
#     copyto!(e, x.e)
#     for i in eachindex(y.e)
#         e[i] = max(e[i], y.e[i])
#     end
#     return PFRational{T}(e)
# end

# """
#     scgd(x::PFRational, y::PFRational)
# here we define `scgd` of rational, which means gcd of both numerator and denominator
# """
# scgd(x::PFRational, y::PFRational) = begin
#     T = promote_type(eltype(x.e), eltype(y.e))
#     e = zeros(T, max(length(x.e), length(y.e)))
#     copyto!(e, x.e)
#     for i in eachindex(y.e)
#         e_min, e_max = minmax(e[i], y.e[i])
#         if e_min > 0
#             e[i] = e_min
#         elseif e_max < 0
#             e[i] = e_max
#         else
#             e[i] = 0
#         end
#     end
#     for i in length(y.e)+1:length(x.e)
#         e[i] = 0
#     end
#     return PFRational{T}(e)
# end

# # assume length(y) <= length(x)
# _copy!(x::PFRational{T}, y::PFRational{T}) where {T<:Integer} = begin
#     for i in eachindex(y.e)
#         x.e[i] = y.e[i]
#     end
#     for i in length(y.e)+1:length(x.e)
#         x.e[i] = 0
#     end
# end

# # assume length(y) <= length(x)
# _mul!(x::PFRational{T}, y::PFRational{T}) where {T<:Integer} = begin
#     for i in eachindex(y.e)
#         x.e[i] += y.e[i]
#     end
# end

# # assume length(y) <= length(x)
# _div!(x::PFRational{T}, y::PFRational{T}) where {T<:Integer} = begin
#     for i in eachindex(y.e)
#         x.e[i] -= y.e[i]
#     end
# end

# # assume length(y) <= length(x)
# _gcd!(x::PFRational{T}, y::PFRational{T}) where {T<:Integer} = begin
#     for i in eachindex(y.e)
#         x.e[i] = min(x.e[i], y.e[i])
#     end
# end

# # assume length(y) <= length(x)
# _lcm!(x::PFRational{T}, y::PFRational{T}) where {T<:Integer} = begin
#     for i in eachindex(y.e)
#         x.e[i] = max(x.e[i], y.e[i])
#     end
# end

# # assume length(y) <= length(x)
# _sgcd!(x::PFRational{T}, y::PFRational{T}) where {T<:Integer} = begin
#     for i in eachindex(y.e)
#         e_min, e_max = minmax(x.e[i], y.e[i])
#         if e_min > 0
#             x.e[i] = e_min
#         elseif e_max < 0
#             x.e[i] = e_max
#         else
#             x.e[i] = 0
#         end
#     end
#     for i in length(y.e)+1:length(x.e)
#         x.e[i] = 0
#     end
# end

# # first n primes
# _prime_table = Int[2, 3, 5, 7]
# # first n primes and their powers
# _prime_power_table = Matrix{BigInt}(undef, 0, 0)
# # positive integers in PFRational form
# _pf_integers = PFRational{Int16}[one(PFRational{Int16})]
# # binomial coefficients in PFRational form
# _pf_bin_data =  PFRational{Int16}[one(PFRational{Int16})]
# # maximum n for binomial coefficients
# _pf_bin_nmax = Int(0)

# function _extend_prime_table(n::Integer)
#     if n <= length(_prime_table)
#         return
#     end
#     start = length(_prime_table)+1
#     for i in start:n
#         push!(_prime_table, nextprime(_prime_table[end] + 1))
#     end
# end

# function _extend_prime_power_table(n::Integer, k::Integer)
#     if n <= size(_prime_power_table, 1) && k <= size(_prime_power_table, 2)
#         return
#     end
#     _extend_prime_table(max(n, size(_prime_power_table, 1)))
#     global _prime_power_table = Matrix{BigInt}(undef, n, k)
#     for i in 1:n
#         for j in 1:k
#             global _prime_power_table[i, j] = big(_prime_table[i])^j
#         end
#     end
# end

# function _extend_pf_integers(n::Integer)
#     if n <= length(_pf_integers)
#         return
#     end
#     while n > _prime_table[end]
#         push!(_prime_table, nextprime(_prime_table[end] + 1))
#     end
#     start = length(_pf_integers)+1
#     for i in start:n
#         f = factor(i)
#         idx = searchsortedfirst(_prime_table, f.pe[end].first)
#         e = zeros(Int16, idx)
#         df = Dict(f.pe)
#         for j in 1:idx
#             e[j] = get(df, _prime_table[j], 0)
#         end
#         push!(_pf_integers, PFRational(e))
#     end
# end


# @inline nth_stored_prime(n::Int)::Int = _prime_table[n]::Int
# @inline get_pf_nmax()::Int = _pf_bin_nmax::Int
# @inline get_pf_bin_data()::Vector{PFRational{Int16}} = _pf_bin_data::Vector{PFRational{Int16}}
# @inline function pf_binomial(n::Integer, k::Integer)
#     0 <= k <= n <= get_pf_nmax() || return one(PFRational{Int16})
#     k = min(k, n-k)
#     return @inbounds get_pf_bin_data()[binomial_index(n, k)]
# end

# @inline function _pf_bin(n::Integer, k::Integer)
#     k = min(k, n-k)
#     return @inbounds get_pf_bin_data()[binomial_index(n, k)]
# end

# function calc_pf_binomial(n::Integer, k::Integer)::PFRational
#     0 <= k <= n || return one(PFRational)
#     k = min(k, n-k)
#     r = one(PFRational{Int16})
#     for i in 1:k
#         r *= _pf_integers[n-i+1] / _pf_integers[i]
#     end
#     return r
# end

# function _extend_pf_binomial(n::Integer)
#     if n <= _pf_bin_nmax
#         return
#     end
#     _extend_pf_integers(n)
#     old_data = copy(_pf_bin_data)
#     resize!(_pf_bin_data, binomial_data_size(n))
#     copyto!(_pf_bin_data, old_data)
#     for m = _pf_bin_nmax+1:n
#         for k = 0:div(m, 2)
#             global _pf_bin_data[binomial_index(m, k)] = calc_pf_binomial(m, k)
#         end
#     end
#     global _pf_bin_nmax = n
# end

# Base.numerator(x::PFRational) = begin
#     num = one(BigInt)
#     for i in eachindex(x.e)
#         if x.e[i] > 0
#             num *= _prime_power_table[i, x.e[i]]
#         end
#     end
#     return num
# end

# Base.denominator(x::PFRational) = begin
#     den = one(BigInt)
#     for i in eachindex(x.e)
#         if x.e[i] < 0
#             den *= _prime_power_table[i, -x.e[i]]
#         end
#     end
#     return den
# end

# Base.convert(::Type{Rational}, x::PFRational) = begin
#     num = one(BigInt)
#     den = one(BigInt)
#     for i in eachindex(x.e)
#         if x.e[i] > 0
#             num *= _prime_power_table[i, x.e[i]]
#         elseif x.e[i] < 0
#             den *= _prime_power_table[i, -x.e[i]]
#         end
#     end
#     return Rational(num, den)
# end

# Base.show(io::IO, x::PFRational) = begin
#     show(io, convert(Rational, x))
# end

# function common_factor!(xs::Vector{PFRational{T}}) where {T<:Integer}
#     if length(xs) == 0
#         return one(PFRational{T})
#     end
#     cf = xs[1]
#     for i in 2:length(xs)
#         cf = gcd(cf, xs[i])
#     end
#     for i in 1:length(xs)
#         xs[i] /= cf
#     end
#     return cf
# end

# function _pf_CG(dj1::Int, dj2::Int, dj3::Int, dm1::Int, dm2::Int, dm3::Int)
#     check_jm(dj1, dm1) && check_jm(dj2, dm2) && check_jm(dj3, dm3) || return zero(SqrtRational)
#     check_couple(dj1, dj2, dj3) || return zero(SqrtRational)
#     dm1 + dm2 == dm3 || return zero(SqrtRational)
#     J::Int = div(dj1 + dj2 + dj3, 2)
#     Jm1::Int = J - dj1
#     Jm2::Int = J - dj2
#     Jm3::Int = J - dj3
#     j1mm1::Int = div(dj1 - dm1, 2)
#     j2mm2::Int = div(dj2 - dm2, 2)
#     j3mm3::Int = div(dj3 - dm3, 2)
#     j2pm2::Int = div(dj2 + dm2, 2)
    
#     Jsize::Int = searchsortedfirst(_prime_table, J + 1)
#     A = pf_alloc(Int16, Jsize)
#     _copy!(A, _pf_bin(dj1, Jm2))
#     _mul!(A, _pf_bin(dj2, Jm3))
#     _div!(A, _pf_bin(J + 1, Jm3))
#     _div!(A, _pf_bin(dj1, j1mm1))
#     _div!(A, _pf_bin(dj2, j2mm2))
#     _div!(A, _pf_bin(dj3, j3mm3))

#     low::Int = max(zero(Int), j1mm1 - Jm2, j2pm2 - Jm1)
#     high::Int = min(Jm3, j1mm1, j2pm2)
#     Bs = Vector{PFRational{Int16}}(undef, high - low + 1)
#     for z in low:high
#         Bs[z-low+1] =  pf_alloc(Int16, Jsize)
#         _copy!(Bs[z-low+1], _pf_bin(Jm3, z))
#         _mul!(Bs[z-low+1], _pf_bin(Jm2, j1mm1 - z))
#         _mul!(Bs[z-low+1], _pf_bin(Jm1, j2pm2 - z))
#     end
#     fB, B = stagger_sum!(Bs)
#     num = one(BigInt)
#     den = one(BigInt)
#     for i in eachindex(A.e)
#         fB.e[i] += div(A.e[i], 2)
#         A.e[i] = rem(A.e[i], 2)
#         if (fB.e[i] > 0 && A.e[i] < 0) || (fB.e[i] < 0 && A.e[i] > 0)
#             A.e[i] = -A.e[i]
#             fB.e[i] = fB.e[i] - A.e[i]
#         end
#         @show i, A.e[i]
#         if A.e[i] > 0
#             num *= _prime_table[i]
#         elseif A.e[i] < 0
#             den *= _prime_table[i]
#         end
#     end
#     rA = num//den
#     rB = convert(Rational, fB) * B
#     return SqrtRational(iphase(low) * rB, rA)
# end

# function stagger_sum!(xs::Vector{PFRational{Int16}})
#     cf = common_factor!(xs)
#     sum = zero(BigInt)
#     sign = 1
#     for i in 1:length(xs)
#         if sign > 0
#             sum += numerator(xs[i])
#         else
#             sum -= numerator(xs[i])
#         end
#         sign = -sign
#     end
#     return cf, sum
# end