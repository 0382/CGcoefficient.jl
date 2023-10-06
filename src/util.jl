# this file contains some useful functions

"""
    iphase(n::Integer)
``(-1)^n``
"""
@inline iphase(n::T) where {T <: Integer} = iseven(n) ? one(T) : -one(T)

"""
    is_same_parity(x::T, y::T) where {T <: Integer}
judge if two integers are same odd or same even
"""
@inline is_same_parity(x::T, y::T) where {T <: Integer} = iseven(x ⊻ y)

"""
    check_jm(dj::T, dm::T) where {T <: Integer}
check if the m-quantum number if one of the components of the j-quantum number,
in other words, `m` and `j` has the same parity, and `abs(m) < j`
"""
@inline check_jm(dj::T, dm::T) where {T <: Integer} = is_same_parity(dj, dm) & (abs(dm) <= dj)

"""
    check_couple(dj1::T, dj2::T, dj3::T) where {T <: Integer}
check if three angular monentum number `j1, j2, j3` can couple
"""
@inline check_couple(dj1::T, dj2::T, dj3::T) where {T <: Integer} = begin
    (dj1 >= 0) & (dj2 >= 0) & is_same_parity(dj1 + dj2, dj3) & (abs(dj1 - dj2) <= dj3 <= dj1 + dj2)
end

@doc raw"""
    binomial_data_size(n::Int)::Int
Store (half) binomial data in the order of
```text
# n - data
  0   1
  1   1
  2   1 2
  3   1 3
  4   1 4 6
  5   1 5 10
  6   1 6 15 20
```
Return the total number of the stored data.
"""
@inline function binomial_data_size(n::Int)::Int
    x = div(n, 2) + 1
    return x * (x + isodd(n))
end

@doc raw"""
    binomial_index(n::Int, k::Int)::Int
Return the index of the binomial coefficient in the table.
"""
@inline function binomial_index(n::Int, k::Int)::Int
    x = div(n, 2) + 1
    return x * (x - iseven(n)) + k + 1
end