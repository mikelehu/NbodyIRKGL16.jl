
struct VecArray{s_,T,dim}
    data::Array{T,dim}
end

@inline function Base.getindex(v::VecArray{s,T,dim},i...) where {s,T,dim}
        Vec{s,T}(NTuple{s,T}(@inbounds v.data[k,i...] for k=1:s))
end

@inline function Base.setindex!(v::VecArray{s,T,dim},vi::Vec{s,T},i...) where {s,T,dim}
    @inbounds for j in 1:s
        v.data[j,i...] = vi[j]
    end
    return nothing
end

@inline function Base.setindex!(v::VecArray{s,T,dim},vi::T2,i...) where {s,T,T2,dim}
    vi_ = convert(T,vi)
    @inbounds for k in 1:s
        v.data[k,i...] = vi_
    end
    return nothing
end

@inline function getindex_(v::VecArray{s,T,dim},i::Int64) where {s,T,dim}
      j = s*(i-1)
      Vec{s,T}(NTuple{s,T}(@inbounds v.data[k+j] for k=1:s))
end


@inline function setindex_!(v::VecArray{s,T,dim},vi::Vec{s,T},i::Int64) where {s,T,dim}
    j = s*(i-1)
    @inbounds for k in 1:s
        v.data[j+k] = vi[k]
    end
    return nothing
end

@inline function setindex_!(v::VecArray{s,T,dim},vi::T2,i::Int64) where {s,T,T2,dim}
    vi_ = convert(T,vi)
    j = s*(i-1)
    @inbounds for k in 1:s
        v.data[j+k] = vi_
    end
    return nothing
end
