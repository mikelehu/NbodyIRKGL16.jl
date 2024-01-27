function InitialPythagorean3BPREL(T=Float64)

    Gm = [parse(BigFloat,"5"), parse(BigFloat,"4"),parse(BigFloat,"3")]

    N=length(Gm)
    x1 = parse(BigFloat,"1")
    y1 = parse(BigFloat,"-1")
#    z1 = parse(BigFloat,"0")
    x2 = parse(BigFloat,"-2")
    y2 = parse(BigFloat,"-1")
#    z2 = parse(BigFloat,"0")
    x3 = parse(BigFloat,"1")
    y3 = parse(BigFloat,"3")
#    z3 = parse(BigFloat,"0")

#    q =  [x1,y1,z1,x2,y2,z2,x3,y3,z3]
    q =  [x2-x1,y2-y1,x3-x2,y3-y2,x1-x3,y1-y3]
    v = zeros(BigFloat,size(q))
    
#    q0 = reshape(q,3,:)
#    v0 = reshape(v,3,:)
    
    q0 = reshape(q,2,:)
    v0 = reshape(v,2,:)
    
    u0 = Array{T}(undef,2,N,2)
    u0[:,:,2] .= convert.(T,v0)
    u0[:,:,1] .= convert.(T,q0)
    bodylist = ["Body-1", "Body-2", "Body-3"]

    return u0, convert.(T,Gm), bodylist

end
