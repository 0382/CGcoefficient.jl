# this file contains some useful functions

"""
    iphase(n::Integer)
``(-1)^n``
"""
iphase(n::T) where {T <: Integer} = iseven(n) ? one(T) : -one(T)

"""
    is_same_parity(x::T, y::T) where {T <: Integer}
judge if two integers are same odd or same even
"""
is_same_parity(x::T, y::T) where {T <: Integer} = iseven(x ⊻ y)

"""
devide integer `n` with 2, and return how to show it,
for odd number, show `"n/2"`, else show the exact value.
"""
show_half(x::Integer) = iseven(x) ? "$(x >> 1)" : "$x/2"

"""
    check_jm(dj::T, dm::T) where {T <: Integer}
check if the m-quantum number if one of the components of the j-quantum number,
in other words, `m` and `j` has the same parity, and `abs(m) < j`
"""
check_jm(dj::T, dm::T) where {T <: Integer} = is_same_parity(dj, dm) & (abs(dm) <= dj)

"""
if `check_jm(j, m)` failed, show this message
"""
jm_mismatching_msg(dj::T, dm::T) where {T <: Integer} = begin
    "j = $(show_half(dj)), m = $(show_half(dm)) are not matched" 
end

"""
    check_couple(dj1::T, dj2::T, dj3::T) where {T <: Integer}
check if three angular monentum number `j1, j2, j3` can couple
"""
check_couple(dj1::T, dj2::T, dj3::T) where {T <: Integer} = begin
    (dj1 >= 0) & (dj2 >= 0) & is_same_parity(dj1 + dj2, dj3) & (abs(dj1 - dj2) <= dj3 <= dj1 + dj2)
end

"""
if `check_couple(j1, j2, j3)` failed, show this message
"""
miss_couple_msg(dj1::T, dj2::T, dj3::T) where {T <: Integer} = begin
    "$(show_half(dj1)), $(show_half(dj2)) cannnot couple to $(show_half(dj3))"
end