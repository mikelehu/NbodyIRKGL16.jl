#   
#  IRKNGL Step functions

#      IRKNGLstep_SIMD_fixed! 


function IRKNGLstep_SIMD_fixed!(ttj, uj,ej, dts,stats, irknglcache::IRKNGL_SIMD_Cache{floatT,fType,pType,s_,dim_}) where {floatT,fType,pType,s_,dim_}


    f = irknglcache.odef
    p = irknglcache.p
    b = irknglcache.b
    c = irknglcache.c
    mu = irknglcache.mu
    nu = irknglcache.nu
    U = irknglcache.U
    U_ = irknglcache.U_
    L = irknglcache.L
    F = irknglcache.F
    Dmin = irknglcache.Dmin
    step_number = irknglcache.step_number[]
    initial_extrap = irknglcache.initial_extrap[]
    len = irknglcache.length_u
    lenq = irknglcache.length_q
    tf= irknglcache.tf
    step_retcode = irknglcache.step_retcode


    s = length(b)
    len = length(uj)
    maxiters = (step_number==1 ? 10+irknglcache.maxiters : irknglcache.maxiters )
    tj = ttj[1]
    te = ttj[2]


    indices=1:len
    indices1 = 1:lenq
    indices2 = (lenq+1):len

    dt=dts[1]
    dtprev=dts[2]
    signdt=dts[3]
    sdt=dt*signdt

    j_iter = 0  # counter of fixed_point iterations
    nf=0
    nf2=0

    if (initial_extrap==1)

      for k in indices2
         Lk = getindex_(L,k)
         dUk = muladd(nu[1], Lk[1], ej[k])
         for is in 2:s
             dUk = muladd(nu[is], Lk[is], dUk)
         end
         setindex_!(U, uj[k]+dUk, k)
      end

    else
    
       for k in indices2
         uej = uj[k] + ej[k]
         setindex_!(U, uej, k)
       end
       
    end

    nf+=s
    f(F, U, p, tj + sdt*c, 1)
    
    for k in indices1
        Fk = getindex_(F,k)
        Lk =sdt*(b*Fk)
        setindex_!(L, Lk, k)
        dUk = muladd(mu[1], Lk[1], ej[k])
        for is in 2:s
            dUk = muladd(mu[is], Lk[is], dUk)
        end
        setindex_!(U, uj[k]+dUk, k)
    end

    Dmin .= Inf

    iter = true # Initialize iter outside the for loop
    plusIt=true
    diffU = true

    @inbounds while (j_iter<maxiters && iter)

         iter = false
         j_iter += 1

         U_.data .= U.data

         nf2+=s
         f(F, U, p, tj + sdt*c, 2)

         for k in indices2
                 Fk = getindex_(F,k)
                 Lk =sdt*(b*Fk)
                 setindex_!(L, Lk, k)
                 dUk = muladd(mu[1], Lk[1], ej[k])
                 for is in 2:s
                     dUk = muladd(mu[is], Lk[is], dUk)
                 end
                 setindex_!(U, uj[k] + dUk, k)
         end

         nf+=s
         f(F, U, p, tj + sdt*c, 1)

         for k in indices1
                 Fk = getindex_(F,k)
                 Lk = sdt*(b*Fk)
                 setindex_!(L, Lk, k)
                 dUk = muladd(mu[1], Lk[1], ej[k])
                 for is in 2:s
                     dUk = muladd(mu[is], Lk[is], dUk)
                 end
                 setindex_!(U, uj[k]+dUk, k)
         end


        diffU = false

        for k in indices1   

             Uk = getindex_(U,k)
             Uk_ = getindex_(U_,k)
             DY = maximum(abs(Uk-Uk_))

             if DY>0
                 diffU = true
                 if DY< Dmin[k]
                    Dmin[k]=DY
                    iter=true
                 end
             end
         end

         if (!iter && diffU && plusIt)  #
             iter=true
             plusIt=false
         else
             plusIt=true
         end

    end # while

    if (!iter && (j_iter==maxiters))
        println("Fail !!!  Step=",tj+te, " dt=", dts[1])
        step_retcode=Failure
        return nothing
    end

    @inbounds if (j_iter<maxiters && diffU)

         j_iter += 1

         nf2+=s
         f(F, U, p, tj + sdt*c, 2)

         for k in indices2
                 Fk = getindex_(F,k)
                 Lk = sdt*(b*Fk)
                 dUk = muladd(mu[1], Lk[1], ej[k])
                 for is in 2:s
                     dUk = muladd(mu[is], Lk[is], dUk)
                 end
                 setindex_!(U, uj[k]+dUk, k)
                 setindex_!(L, Lk, k)
         end

         nf+=s
         f(F, U, p, tj + sdt*c, 1)

         for k in indices1
                 Fk = getindex_(F,k)
                 Lk = sdt*(b*Fk)
                 setindex_!(L, Lk, k)
         end

    end


    @inbounds for k in indices    #Equivalent to compensated summation
         Lk = getindex_(L,k)
         L_sum = sum(Lk)
         res = Base.TwicePrecision(uj[k], ej[k]) + L_sum
         uj[k] = res.hi
         ej[k] = res.lo
    end

    res = Base.TwicePrecision(tj, te) + sdt
    ttj[1] = res.hi
    ttj[2] = res.lo

    dts[1]=min(abs(dt),abs(tf-(ttj[1]+ttj[2])))
    dts[2]=dt

    if dts[2]!=Inf
        stats.nnonliniter+=j_iter
        stats.nf+=nf
        stats.nf2+=nf2
    end

    retcode=Success
    return nothing

 end


