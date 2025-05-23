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
@inline is_same_parity(x::T, y::T) where {T <: Integer} = mod(x, 2) == mod(y, 2)

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
    iseven(dj1 + dj2 + dj3) && (dj1 <= dj2 + dj3) && (dj2 <= dj1 + dj3) && (dj3 <= dj1 + dj2)
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

"""
    check_CG(dj1::T, dj2::T, dj3::T, dm1::T, dm2::T, dm3::T) where {T <: Integer}
Check if the Clebsch-Gordan coefficient is valid.
"""
@inline function check_CG(dj1::T, dj2::T, dj3::T, dm1::T, dm2::T, dm3::T) where {T <: Integer}
    check_jm(dj1, dm1) && check_jm(dj2, dm2) && check_jm(dj3, dm3) && check_couple(dj1, dj2, dj3) && (dm1 + dm2 == dm3)
end

"""
    check_Gaunt(l1::T, l2::T, l3::T, m1::T, m2::T, m3::T) where {T <: Integer}
Check if the Gaunt coefficient is valid.
"""
@inline function check_Gaunt(l1::T, l2::T, l3::T, m1::T, m2::T, m3::T) where {T <: Integer}
    if !check_couple(l1, l2, l3)
        return false
    end
    if abs(m1) > l1 || abs(m2) > l2 || abs(m3) > l3
        return false
    end
    if m1 + m2 + m3 != 0
        return false
    end
    return true
end

"""
    check_3j(dj1::T, dj2::T, dj3::T, dm1::T, dm2::T, dm3::T) where {T <: Integer}
Check if the Wigner 3-j symbol is valid.
"""
@inline function check_3j(dj1::T, dj2::T, dj3::T, dm1::T, dm2::T, dm3::T) where {T <: Integer}
    check_jm(dj1, dm1) && check_jm(dj2, dm2) && check_jm(dj3, dm3) && check_couple(dj1, dj2, dj3) && (dm1 + dm2 + dm3 == 0)
end

"""
    check_6j(dj1::T, dj2::T, dj3::T, dj4::T, dj5::T, dj6::T) where {T <: Integer}
Check if the Wigner 6-j symbol is valid.
"""
@inline function check_6j(dj1::T, dj2::T, dj3::T, dj4::T, dj5::T, dj6::T) where {T <: Integer}
    check_couple(dj1, dj2, dj3) && check_couple(dj1, dj5, dj6) && check_couple(dj4, dj2, dj6) && check_couple(dj4, dj5, dj3)
end

"""
    check_9j(dj1::T, dj2::T, dj3::T, dj4::T, dj5::T, dj6::T, dj7::T, dj8::T, dj9::T) where {T <: Integer}
Check if the Wigner 9-j symbol is valid.
"""
@inline function check_9j(dj1::T, dj2::T, dj3::T, dj4::T, dj5::T, dj6::T, dj7::T, dj8::T, dj9::T) where {T <: Integer}
    check_couple(dj1, dj2, dj3) && check_couple(dj4, dj5, dj6) && check_couple(dj7, dj8, dj9) &&
    check_couple(dj1, dj4, dj7) && check_couple(dj2, dj5, dj8) && check_couple(dj3, dj6, dj9)
end

"""
    check_Moshinsky(N::T, L::T, n::T, l::T, n1::T, l1::T, n2::T, l2::T, Λ::T) where {T <: Integer}
Check if the Moshinsky coefficient is valid.
"""
@inline function check_Moshinsky(N::T, L::T, n::T, l::T, n1::T, l1::T, n2::T, l2::T, Λ::T) where {T <: Integer}
    if N < 0 || n < 0 || n1 < 0 || n2 < 0
        return false
    end
    if Λ > L + l || L > Λ + l || l > Λ + L
        return false
    end
    if Λ > l1 + l2 || l1 > Λ + l2 || l2 > Λ + l1
        return false
    end
    E = 2N + L
    e = 2n + l
    e1 = 2n1 + l1
    e2 = 2n2 + l2
    if E + e != e1 + e2
        return false
    end
    return true
end