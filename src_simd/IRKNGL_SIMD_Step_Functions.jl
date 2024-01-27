#   
#  IRKNGL Step functions

#      IRKNGLstep_SIMD_fixed!
#      IRKNGLstep_SIMD_adap!   


 function IRKNGLstep_SIMD_fixed!(ttj,tj1, uj,ej,prob,dts, irknglcache::IRKNGL_SIMD_Cache{floatType,fType,pType,s_,dim,dim_}) where {floatType,fType,pType,s_,dim,dim_}

    f = irknglcache.odef
    p = irknglcache.p
    b = irknglcache.b
    c = irknglcache.c
    mu = irknglcache.mu
    nu = irknglcache.nu
    U = irknglcache.U
    U_ = irknglcache.U_
    L = irknglcache.L
    dU = irknglcache.dU
    F = irknglcache.F
    Dmin = irknglcache.Dmin
    step_number = irknglcache.step_number[]
    initial_interp = irknglcache.initial_interp[]
    len = irknglcache.length_u
    lenq = irknglcache.length_q
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


    if (initial_interp==1)

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

    j_iter = 0  # counter of fixed_point iterations
    Dmin .= Inf

    iter = true # Initialize iter outside the for loop
    plusIt=true
    diffU = true

    @inbounds while (j_iter<maxiters && iter)

         iter = false
         j_iter += 1

         U_.data .= U.data

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

        for k in indices1   # Hemen indices1 jarri liteke, q'=v, v'=f(q,t) moduko ED-a dela suposatuz

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
        return("Failure",0)
    end

    @inbounds if (j_iter<maxiters && diffU)

         j_iter += 1

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

         f(F, U, p, tj + sdt*c, 1)

         for k in indices1
                 Fk = getindex_(F,k)
                 Lk = sdt*(b*Fk)
                 setindex_!(L, Lk, k)
         end

    end


     @inbounds for k in indices    #Batura konpentsatuaren parekoa
         Lk = getindex_(L,k)
         L_sum = sum(Lk)
         res = Base.TwicePrecision(uj[k], ej[k]) + L_sum
         uj[k] = res.hi
         ej[k] = res.lo
      end

      res = Base.TwicePrecision(tj, te) + sdt
      ttj[1] = res.hi
      ttj[2] = res.lo

      dts[1]=min(abs(dt),abs(tj1-(ttj[1]+ttj[2])))
      dts[2]=dt


      return  ("Success",j_iter)

 end




function IRKNGLstep_SIMD_adap!(ttj,tj1, uj,ej,prob,dts,irkngl_cache::IRKNGL_SIMD_Cache{floatType,fType,pType,s_,dim,dim_}) where {floatType,fType,pType,s_,dim,dim_}

    trace = false

    f = irkngl_cache.odef
    p = irkngl_cache.p
    b = irkngl_cache.b
    c = irkngl_cache.c
    mu = irkngl_cache.mu
    nu = irkngl_cache.nu
    U = irkngl_cache.U
    U_ = irkngl_cache.U_
    L = irkngl_cache.L
    dU = irkngl_cache.dU
    theta=irkngl_cache.theta
    omega=irkngl_cache.omega
    alpha=irkngl_cache.alpha
    g=irkngl_cache.g
    d=irkngl_cache.d
    s_beta=irkngl_cache.s_beta
    F = irkngl_cache.F
    Dmin = irkngl_cache.Dmin
    step_number = irkngl_cache.step_number[]
    tau= irkngl_cache.tau[]
    initial_interp = irkngl_cache.initial_interp[]
    len = irkngl_cache.length_u
    lenq = irkngl_cache.length_q
    nrmdigits=irkngl_cache.nrmdigits[]
    E_weights=irkngl_cache.E_weights

    s = length(b)
    maxiters = (step_number==1 ? 10+irkngl_cache.maxiters : irkngl_cache.maxiters )
    tj = ttj[1]
    te = ttj[2]


    indices=1:len
    indices1 = 1:lenq
    indices2 = (lenq+1):len

    uEtype=eltype(uj)
    R=uEtype(2)^(precision(uEtype)-16)

    dt=dts[1]
    dt_=dts[2]
    signdt=dts[3]  
    sdt=dt*signdt

    if (initial_interp==1)

        for k in indices2
            Lk = getindex_(L,k)
            dUk = muladd(nu[1], Lk[1], ej[k])
            for js in 2:s
                dUk = muladd(nu[js], Lk[js], dUk)
            end
            setindex_!(dU, dUk, k)
        end
        
        lambda=dt/dt_-1

        if abs(lambda)>eps(uEtype) 
            

            for j in 1:s+1

                thetaj=getindex_(theta,j)
                omegaj=getindex_(omega,j)
                betaj=1/muladd(lambda,thetaj,omegaj)
                setindex_!(alpha, betaj,j)

            end

            for i in 1:s_
                s_betai=1/sum(alpha.data[i,1:s+1])
                s_beta=Base.setindex(s_beta,s_betai,i)  # s_beta[i]=s_betai
#                if isinf(s_betai)
#                    println("s_betai infinitua !!!!!!!!!!!!!!!!!!!!!!")
#                end
            end

            for j in 1:s+1
               alphaj=getindex_(alpha,j)
               setindex_!(alpha, alphaj*s_beta, j)
            end
  
            for k in indices2
                    dUk=getindex_(dU,k)
                    dUik = alpha[1]*ej[k]
                    for js in 1:s
                        dUik = muladd(alpha[js+1], dUk[js], dUik)
                    end
                    ujdU=uj[k]+dUik
                    setindex_!(U, ujdU, k)
            end              
            
        else    

            for k in indices2
                dUk=getindex_(dU,k)
                ujdU=uj[k]+dUk
                setindex_!(U, ujdU, k)
            end
            
        end  

    else

        for k in indices2
            uej = uj[k] + ej[k]
           setindex_!(U, uej, k)
         end

    end
    
    f(F, U, p, tj + sdt*c, 1 )

    for k in indices1
        Fk = getindex_(F,k)
        Lk = sdt*(b*Fk)
        dUk = muladd(mu[1], Lk[1], ej[k])
        for js in 2:s
            dUk = muladd(mu[js], Lk[js], dUk)
        end
        Uk = uj[k]+dUk
        Uk_=getindex_(U,k)
        setindex_!(U_, Uk_,k)
        setindex_!(U, Uk, k)
        setindex_!(L, Lk, k)    
    end
    
    trace ? println("Urratsa tj=",tj+te,", iter=", 1, ",dt=", Float32(dt)) : nothing

    j_iter = 0
    Dmin .= Inf
    iter = true # Initialize iter outside the for loop
    plusIt=true
    diffU =true
    dtj1=dt
    dtaux=dt
    fase3=false
    Dtmin=Inf
    plusDt=2

    DE=zero(uEtype)
    E2=zero(uEtype)

   @inbounds while (iter && j_iter<maxiters) 

        j_iter += 1
        dtj0=dtj1

        trace ? println("*** ", j_iter, ". iterazioa:") : nothing

        f(F, U, p, tj + sdt*c, 2 )

        for k in indices2
            Fk = getindex_(F,k)
            Lk = sdt*(b*Fk)
            dUk = muladd(mu[1], Lk[1], ej[k])
            for js in 2:s
                dUk = muladd(mu[js], Lk[js], dUk)
            end
            Uk = uj[k]+dUk
            setindex_!(U, Uk, k)
            setindex_!(L, Lk, k)
            setindex_!(dU,dUk,k)
        end

        f(F, U, p, tj + sdt*c, 1 )

        diffU = false
        iter= false

        for k in indices1

            Fk = getindex_(F,k)
            Lk = sdt*(b*Fk)
            dUk = muladd(mu[1], Lk[1], ej[k])
            for js in 2:s
                dUk = muladd(mu[js], Lk[js], dUk)
            end
            Uk = uj[k]+dUk
            Uk_=getindex_(U,k)
            setindex_!(U_, Uk_,k)
            setindex_!(U, Uk, k)
            setindex_!(L, Lk, k)

            # DY = maximum(abs(Uk-Uk_))
            diffk=Rdigits(Uk,nrmdigits)-Rdigits(Uk_,nrmdigits)
            DY = maximum(abs(diffk))

            if DY>0
                diffU = true
                if DY< Dmin[k]
                   Dmin[k]=DY
                   iter=true
                end
            end
        
        end

        if (!iter && diffU && plusIt)
            iter=true
            plusIt=false
        else
            plusIt=true
        end

        if trace == true
            Dminaux=replace(Dmin, Inf=>0)
            println("iter=", iter, ",diffU=", diffU, ",plusIt=", plusIt, ", norm(Dmin)=",norm(Dminaux)) 
        end

        E2=zero(uEtype)
        DE=zero(uEtype)

        for k in indices2 

            Fk = getindex_(F,k)

            Ek=dt*sum(g*Fk)   # Ek =dt*dot(g, Fk)
            Dk=sum(d*Fk)      # Dk=  dot(d, Fk) 
            wEk= E_weights[k]*Ek
            E2=muladd(Ek,wEk,E2)
            DE=muladd(Dk,wEk,DE)

        end

        fz = log(E2) - (2s-2)*log(tau)
        fderz = (2s-2)+dt*DE/E2
        dtaux = Rdigits(dt*exp(-fz/fderz),R)
        
        dtj1=min(dtaux,abs(tj1-(ttj[1]+ttj[2])))
        Dt=abs((dtj0-dtj1)/dtj1)

        if ((Dtmin<1e-6) || !iter) fase3=true end
            
        if (0<Dt<1e-1) && !fase3

            if (Dt>Dtmin && plusDt==1) 
                fase3=true 
                plusDt-=1
            elseif (Dt>Dtmin && plusDt>1) 
                  plusDt-=1 
                else 
                  plusDt=2
            end  


            if trace
                    println("Dt<0.0001", " fase3=",fase3, ", dt=",Float32(dt) ,", dtj0=", Float32(dtj0),
                " ,dtj1=", Float32(dtj1), ", Dtmin=", Float32(Dtmin),", Dt=", Float32(Dt), ", plusDt=", plusDt)
                    
            end


            if Dt<Dtmin 

               if trace 
                  println("eguneratu", ",dt=",Float32(dtj1) )
                  println()
               end

               Dtmin=Dt
          
               lambda=dtj1/dt-1

               #= Nola bektorizatu?=#

               for j in 1:s+1

                thetaj=getindex_(theta,j)
                omegaj=getindex_(omega,j)
                betaj=1/muladd(lambda,thetaj,omegaj)
                setindex_!(alpha, betaj,j)

               end

               for i in 1:s_
                  s_betai=1/sum(alpha.data[i,1:s+1])
                  s_beta=Base.setindex(s_beta,s_betai,i)  # s_beta[i]=s_betai
               end

               for j in 1:s+1
                 alphaj=getindex_(alpha,j)
                 setindex_!(alpha, alphaj*s_beta, j)
               end


                for k in indices2

                    dUk=getindex_(dU,k)
                    dUik = alpha[1]*ej[k]
                    for js in 1:s
                       dUik = muladd(alpha[js+1], dUk[js], dUik)
                    end
                    ujdU=uj[k]+dUik
                    setindex_!(U, ujdU, k)

                end
           
                dt=dtj1
                sdt=dt*signdt

                f(F, U, p, tj + sdt*c, 1 )

                for k in indices1

                    Fk = getindex_(F,k)
                    Lk = sdt*(b*Fk)
                    dUk = muladd(mu[1], Lk[1], ej[k])
                    for js in 2:s
                        dUk = muladd(mu[js], Lk[js], dUk)
                    end
                    Uk = uj[k]+dUk
                    Uk_=getindex_(U,k)
                    setindex_!(U_, Uk_,k)
                    setindex_!(U, Uk, k)
                    setindex_!(L, Lk, k)
                
                end    

            else
                if trace 
                    println("mantendu", ",dt=",Float32(dt) )
                    println()
                 end
            end

        else  
            if !fase3 
                plusDt=2 
                Dtmin=min(Dt,Dtmin)
            if trace         
                println("Dt>0.0001", " fase3=",fase3, ", dt=",Float32(dt) ,", dtj0=", Float32(dtj0)," ,dtj1=", Float32(dtj1),
                ", Dtmin=", Float32(Dtmin),", Dt=", Float32(Dt), ", plusDt=", plusDt)
                println("mantendu", ",dt=",Float32(dt) )
                println()
            end
            else
                if trace
                println("fase3=", fase3, ",dt=",Float32(dt), " finkoarekin")
                println()
                end
            end

        end


  end   # while
  
 
  if (iter && j_iter==maxiters) 
      println("Failure !!!. Maximum number of iterations.  Step=",tj+te, " dt=", dt)
      return("Failure",0)
  end


  @inbounds if (j_iter<maxiters && diffU)   
      j_iter += 1

      f(F, U, p, tj + sdt*c, 2)

      for k in indices2
          Fk = getindex_(F, k)
          Lk = sdt*(b*Fk)

          dUk = muladd(mu[1], Lk[1], ej[k])
          for js in 2:s
              dUk = muladd(mu[js], Lk[js], dUk)
          end
          Uk = uj[k]+dUk
          setindex_!(U, Uk, k)
          setindex_!(L, Lk, k)
      end

      f(F, U, p, tj + sdt*c, 1 )

      for k in indices1
          Fk = getindex_(F,k)
          Lk = sdt*(b*Fk)
          setindex_!(L, Lk, k)
      end

   end 


  @inbounds for k in indices    #Batura konpentsatuaren parekoa

    Lk = getindex_(L,k)
    L_sum = sum(Lk)

#    res = Base.TwicePrecision(uj[k], ej[k]) + L_sum
    res = Base.TwicePrecision(uj[k], ej[k]) + Rdigits(L_sum,2*nrmdigits)
    uj[k] = res.hi
    ej[k] = res.lo
  end

  res = Base.TwicePrecision(tj, te) + sdt
  ttj[1] = res.hi
  ttj[2] = res.lo

  dts[2]=dt

  # hurrengo urratserako abiapuntua

  dH=-DE/((s-1)*E2)*dt
  trace ?  println("dH=", dH, ",DE=",DE, ",E2=",E2, ",dt=",dt)  : nothing
  if abs(dH)<2
    dt=(2+dH)/(2-dH)*dt
    trace ?  println("Abiapuntu berria h^{0}=", dt)  : nothing 
  end

#  dts[1]=sdt*min(abs(dt),abs(tf-(ttj[1]+ttj[2])))
  dts[1]=min(abs(dt),abs(tj1-(ttj[1]+ttj[2])))

#=
      println("Urratsa: dt_z=",dts[2], " dt_b=",dts[1])
       println("urratsa tj=",tj+te)
       println("")
       println((step_number,Float64(-1+E2/tau^(2s-2))))
    =#
  return  ("Success",j_iter)


end

