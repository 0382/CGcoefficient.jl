mutable struct PrimePowItem
    size::Int
    offset::Int
    ulongdata::Vector{Culong}
    bigdata::Vector{BigInt}
    function PrimePowItem(size::Int, offset::Int, ulongdata::Vector{Culong}, bigdata::Vector{BigInt})
        new(size, offset, ulongdata, bigdata)
    end
end

function make_pow_item(p::Culong, k::Int)
    offset = floor(Int, log(big(p), big(typemax(Culong))))
    ulongdata = Vector{Culong}(undef, offset)
    for i in 1:offset
        ulongdata[i] = p^i
    end
    if k > offset
        bigdata = [BigInt(p) for i in offset+1:k]
        for i in offset+1:k
            MPZ.pow_ui!(bigdata[i - offset], i)
        end
    else
        bigdata = BigInt[]
    end
    PrimePowItem(max(k, offset), offset, ulongdata, bigdata)
end

# first n primes
_prime_table = Culong[2, 3, 5, 7]
# first n primes and their powers
_prime_pow_table = PrimePowItem[]

@inline nth_stored_prime(n::Integer)::Culong = @inbounds _prime_table[n]::Culong
@inline last_stored_prime()::Culong = _prime_table[end]::Culong
@inline stored_prime_num()::Int = length(_prime_table::Vector{Culong})
@inline prime_pow_table()::Vector{PrimePowItem} = _prime_pow_table::Vector{PrimePowItem}
@inline prime_pow_table_size()::Int = length(prime_pow_table())
@inline max_prime_index(n::Integer)::Int = searchsortedlast(_prime_table::Vector{Culong}, convert(Culong, n))

function _extend_prime_table(size::Integer)
    size <= stored_prime_num() && return
    start = stored_prime_num() + 1
    p = last_stored_prime()
    for i in start:size
        p = nextprime(p + 1)
        push!(_prime_table, p)
    end
end

function _extend_prime_table_to(n::Integer)
    p = last_stored_prime()
    n <= p && return
    while p < n
        p = nextprime(p + 1)
        push!(_prime_table, p)
    end
end

function _extend_item!(item::PrimePowItem, k::Integer)
    k <= item.size && return
    k <= item.offset && return
    resize!(item.bigdata, k - item.offset)
    for i in item.size+1:k
        p = item.ulongdata[1]
        idx = i - item.offset
        item.bigdata[idx] = BigInt(p)
        MPZ.pow_ui!(item.bigdata[idx], i)
    end
    item.size = k
end


function _extend_prime_pow_table(size::Integer, maxnum::Integer, times::Int)
    size <= prime_pow_table_size() && return
    _extend_prime_table(size)
    factor = log2(maxnum) * times
    for i in eachindex(prime_pow_table())
        p = nth_stored_prime(i)
        k = floor(Int, factor / log2(p))
        _extend_item!(prime_pow_table()[i], k)
    end
    sizehint!(prime_pow_table(), size)
    start = prime_pow_table_size() + 1
    for i in start:size
        p = nth_stored_prime(i)
        k = floor(Int, factor / log2(p))
        push!(prime_pow_table(), make_pow_item(p, k))
    end
end

function mul_p_pow_k!(x::BigInt, pidx::Int, k::Int)
    item = prime_pow_table()[pidx]::PrimePowItem
    if k <= item.offset
        return MPZ.mul_ui!(x, @inbounds item.ulongdata[k])
    end
    return MPZ.mul!(x, item.bigdata[k - item.offset])
end
