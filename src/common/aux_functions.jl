
#
#  PolInterp, PolInterp!
#  Rdigits


function PolInterp(X::AbstractVector, Y::AbstractVector, n::Int64, z::Number)

    xtype=eltype(X)

    pz = zero(xtype)

    for i=1:n
        lag = one(xtype)
        for j=1:n
            if (j!=i) 
                lag *= X[i]-X[j]
            end
        end
        liz = 1/lag
        for j=1:n
            if (j!=i) 
                liz *= z-X[j]
            end
        end
        pz += Y[i]*liz
    end
    return pz
end




function PolInterp!(P::AbstractVector, Liz::AbstractVector, X::AbstractVector, Y::AbstractVector, n::Int64, Z::AbstractVector)

    xtype=eltype(X)

    P .= zero(xtype)

    for i=1:n
        lag = one(xtype)
        for j=1:n
            if (j!=i) 
                lag *= X[i]-X[j]
            end
        end
        @. Liz = 1/lag
        for j=1:n
            if (j!=i) 
                @. Liz *= Z-X[j]
            end
        end
        @. P += Y[i]*Liz
    end
    return nothing

end


function PolInterp(X::AbstractVector, Y::AbstractVector, n0::Int64, n1::Int64, z::Number)

    xtype=eltype(X)

    pz = zero(xtype)

    for i=n0:n1

        lag = one(xtype)
        for j=n0:n1
            if (j!=i) 
                lag *= X[i]-X[j]
            end
        end
        liz = 1/lag
        for j=n0:n1
            if (j!=i) 
                liz *= z-X[j]
            end
        end
        pz += Y[i]*liz
    end
    return pz
end


#
# 
# 2022-04-26 Bertsio berria sortu dut SIMD inplementazioarentzat
# behin behinekoa
#

function PolInterp(X, Y, n::Int64, z::Number)

    xtype=eltype(X)

    pz = zero(xtype)
    for i=1:n
        lag = one(xtype)
        for j=1:n
            if (j!=i) 
                lag *= X[i]-X[j]
            end
        end
        liz = 1/lag
        for j=1:n
            if (j!=i) 
                liz *= z-X[j]
            end
        end
        pz += Y[i]*liz
    end
    return pz
end



function Rdigits(x::Real,r::Real)

    mx=r*x
    mxx=mx+x
    return mxx-mx

end