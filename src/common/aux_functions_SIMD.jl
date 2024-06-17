
#
#  PolInterp, PolInterp!
#  Rdigits



#
# 2024-05-23
#

function PolInterp!(P::Vec, Liz::Vec, X::Vec, Y::Vec, n::Int64, Z::Vec)

    xtype=eltype(X)

#    P.data[1:n]=zero(xtype)

    for i in 1:n 
        P=Base.setindex(P, 0, i)
    end

    for i=1:n
        lag = one(xtype)
        for j=1:n
            if (j!=i) 
                lag *= X[i]-X[j]
            end
        end

        for i in 1:n 
            Liz=Base.setindex(Liz, 1/lag, i)
        end        

        for j=1:n
            if (j!=i) 
                Liz *= Z-X[j]
            end
        end

        P += Y[i]*Liz
    end

#    return nothing
    return P

end


