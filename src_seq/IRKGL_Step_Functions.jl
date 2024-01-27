#
#  IRKGLstep_fixed!
#  IRKGLstep_adap!


function IRKGLstep_fixed!(ttj,tj1, uj,ej,prob,dts,irkgl_cache::IRKGL_Cache{uType,tType2,fType,pType}) where {uType,tType2,fType,pType}

    trace = false
    
     f = irkgl_cache.odef
     p = irkgl_cache.p
     b = irkgl_cache.b
     c = irkgl_cache.c
     mu = irkgl_cache.mu
     nu = irkgl_cache.nu
     U = irkgl_cache.U
     U_ = irkgl_cache.U_
     L = irkgl_cache.L
     dU = irkgl_cache.dU
     F = irkgl_cache.F
     Dmin = irkgl_cache.Dmin
     step_number = irkgl_cache.step_number[]
     initial_interp = irkgl_cache.initial_interp[]
     len = irkgl_cache.length_u
     s = length(b)
     maxiters = (step_number==1 ? 10+irkgl_cache.maxiters : irkgl_cache.maxiters )
     tj = ttj[1]
     te = ttj[2]
     indices=1:len
    
     dt=dts[1]
     dtprev=dts[2]
     signdt=dts[3]
     sdt=dt*signdt
    

     if (initial_interp==1)
         for is in 1:s
             for k in indices
                 dUik = muladd(nu[is,1], L[1][k], ej[k])
                 for js in 2:s
                    dUik = muladd(nu[is,js], L[js][k], dUik)
                 end
                 U[is][k] =  uj[k]  + dUik
             end
         end
     else
         for is in 1:s
             for k in indices
                U[is][k] = uj[k] + ej[k]
             end
         end
     end
    
    
    j_iter = 0  # counter of fixed_point iterations
    Dmin .= Inf

    iter = true # Initialize iter outside the for loop
    plusIt=true
    diffU = true
   
    
    trace ? println("Urratsa tj=",tj+te,", iter=", 1, ",dt=", Float32(dt)) : nothing

    
    @inbounds while (j_iter<maxiters && iter)  
    
          iter = false
          j_iter += 1
    
          for is in 1:s
              f(F[is], U[is], p,  tj  + sdt*c[is])
              for k in indices
                  U_[is][k] = U[is][k]
                  L[is][k] = sdt*(b[is]*F[is][k])
              end
          end
          
    
          for is in 1:s
              for k in indices
                 dUik = muladd(mu[is,1], L[1][k], ej[k])
                 for js in 2:s
                    dUik = muladd(mu[is,js], L[js][k], dUik)
                 end
                 U[is][k] =  uj[k] + dUik
              end
          end
    
          diffU = false
    
          for k in indices
              DY = abs(U[1][k]-U_[1][k])
              for is in 2:s 
                  DY=max(abs(U[is][k]-U_[is][k]),DY)
              end 
    
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
    
      end # while
    
      if (!iter && (j_iter==maxiters))
         println("Failure !!!  Step=",tj+te, " dt=", dts[1])
         return("Failure",0)
     end
    
      @inbounds if (j_iter<maxiters && diffU)   
        
              j_iter += 1
    
              for is in 1:s
                  f(F[is], U[is], p,  tj  + sdt*c[is])
                  for k in indices
                      L[is][k] = sdt*(b[is]*F[is][k])
                  end
             end
      end
    
    
      @inbounds for k in indices    #Batura konpentsatuaren parekoa
    
          L_sum = L[1][k]
          for is in 2:s
              L_sum+=L[is][k]
          end
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
    


function IRKGLstep_adap!(ttj,tj1, uj,ej,prob,dts,irkgl_cache::IRKGL_Cache{uType,tType2,fType,pType}) where {uType,tType2,fType,pType}
   
    trace = false

    f = irkgl_cache.odef
    p = irkgl_cache.p
    b = irkgl_cache.b
    c = irkgl_cache.c
    mu = irkgl_cache.mu
    nu = irkgl_cache.nu
    U = irkgl_cache.U
    U_ = irkgl_cache.U_
    L = irkgl_cache.L
    dU = irkgl_cache.dU
    theta=irkgl_cache.theta
    omega=irkgl_cache.omega
    alpha=irkgl_cache.alpha
    g=irkgl_cache.g
    d=irkgl_cache.d
    F = irkgl_cache.F
    Dmin = irkgl_cache.Dmin
    step_number = irkgl_cache.step_number[]
    tau= irkgl_cache.tau[]
    initial_interp = irkgl_cache.initial_interp[]
    len = irkgl_cache.length_u
    nrmdigits=irkgl_cache.nrmdigits[]
    E_weights=irkgl_cache.E_weights

    s = length(b)
    maxiters = (step_number==1 ? 10+irkgl_cache.maxiters : irkgl_cache.maxiters )
    tj = ttj[1]
    te = ttj[2]
    indices=1:len


    uEtype=eltype(uj)
    R=uEtype(2)^(precision(uEtype)-16)   #16

    dt=dts[1]
    dt_=dts[2]
    signdt=dts[3]  
    sdt=dt*signdt

    if (initial_interp==1)

     for is in 1:s
         for k in indices
             dUik = muladd(nu[is,1], L[1][k], ej[k])
             for js in 2:s
                dUik = muladd(nu[is,js], L[js][k], dUik)
             end
             dU[is][k] =  dUik
         end
     end
        
     lambda=dt/dt_-1

      if abs(lambda)>eps(uEtype) 
            
               for i in 1:s
                   sumbetai=zero(uEtype)
                   for j in 1:s+1
                      aux=muladd(lambda,theta[i,j],omega[i,j])
                      alpha[i,j]=1/aux
                     sumbetai+=alpha[i,j]
                   end
                %
                   for j in 1:s+1
                    alpha[i,j]=alpha[i,j]/sumbetai
                   end
               end

               for is in 1:s
                 for k in indices
                    dUik = alpha[is,1]*ej[k]
                    for js in 1:s
                        dUik = muladd(alpha[is,js+1], dU[js][k], dUik)
                    end
                    U[is][k] = uj[k] + dUik
                 end
               end
                
            
        else    
            for is in 1:s
              for k in indices
                 U[is][k] = uj[k] + dU[is][k]
              end
            end 
            
      end  

    else
        for is in 1:s
            for k in indices
               U[is][k] = uj[k] + ej[k]
            end
        end
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

        for is in 1:s
            f(F[is], U[is], p,  tj  + sdt*c[is])
            for k in indices
                L[is][k] = sdt*(b[is]*F[is][k])
            end
        end

        for is in 1:s
            for k in indices
               dUik = muladd(mu[is,1], L[1][k], ej[k])
               for js in 2:s
                  dUik = muladd(mu[is,js], L[js][k], dUik)
               end
               U_[is][k] = U[is][k]
               dU[is][k] = dUik
               U[is][k] =  uj[k] + dUik
            end
        end

        diffU = false
        iter= false 
  
        for k in indices 
#            DY = abs(U[1][k]-U_[1][k])
            diff1k = Rdigits(U[1][k],nrmdigits)-Rdigits(U_[1][k],nrmdigits)
            DY = abs(diff1k)
            for is in 2:s 
#                  DY=max(abs(U[is][k]-U_[is][k]),DY)
                diffisk = Rdigits(U[is][k],nrmdigits)-Rdigits(U_[is][k],nrmdigits)
                DY=max(abs(diffisk),DY)
            end 

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
        for k in indices
            Ek=g[1]*F[1][k]
            Dk=d[1]*F[1][k]
            for js in 2:s
                Ek= muladd(g[js], F[js][k],  Ek)
                Dk= muladd(d[js], F[js][k],  Dk)
            end
            Ek=dt*Ek
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

               for i in 1:s
                   sumbetai=zero(uEtype) 
                   for j in 1:s+1
                      aux=muladd(lambda,theta[i,j],omega[i,j])
                      alpha[i,j]=1/aux
                     sumbetai+=alpha[i,j]
                   end
                %
                   for j in 1:s+1
                    alpha[i,j]=alpha[i,j]/sumbetai
                   end
               end

               for is in 1:s
                 for k in indices
                    dUik = alpha[is,1]*ej[k]
                    for js in 1:s
                        dUik = muladd(alpha[is,js+1], dU[js][k], dUik)
                    end
                    U[is][k] = uj[k] + dUik
                 end
               end
                
             dt=dtj1
             sdt=dtj1*signdt
 
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

  #  println("amaitu da:", "j_iter=", j_iter,",iter=", iter,",diffU=", diffU, ",plusIt=", plusIt)
  
 
  if (iter && j_iter==maxiters) 
      println("Failure !!!. Maximum number of iterations.  Step=",tj+te, " dt=", dt)
      return("Failure",0)
  end


  @inbounds if (j_iter<maxiters && diffU)  
      j_iter += 1
      for is in 1:s
        f(F[is], U[is], p,  tj  + sdt*c[is])
        for k in indices
            L[is][k] = sdt*(b[is]*F[is][k])
        end
      end

  end 


  @inbounds for k in indices    #Batura konpentsatuaren parekoa

    L_sum = L[1][k]
    for is in 2:s
        L_sum+=L[is][k]
    end
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
  if abs(dH)<2
    dt=(2+dH)/(2-dH)*dt
    trace ?  println("Abiapuntu berria h^{0}=", dt)  : nothing  #println("Abiapuntu aurreko urratsa")
  end

  dts[1]=min(abs(dt),abs(tj1-(ttj[1]+ttj[2])))

  return  ("Success",j_iter)


end


