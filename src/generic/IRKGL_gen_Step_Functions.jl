#
#  IRKGLstep_fixed!
#  IRKGLstep_adap!


function IRKGLstep_fixed!(ttj,uj,ej,dts,stats,irkgl_cache::NbodyIRKGL16.IRKGL_Cache{uT,tT,fT,pT}) where {uT,tT,fT,pT}

    
     f = irkgl_cache.odef
     p = irkgl_cache.p
     b = irkgl_cache.b
     c = irkgl_cache.c
     a = irkgl_cache.a
     mu = irkgl_cache.mu
     nu = irkgl_cache.nu
     U = irkgl_cache.U
     U_ = irkgl_cache.U_
     L = irkgl_cache.L
     F = irkgl_cache.F
     Dmin = irkgl_cache.Dmin
     step_number = irkgl_cache.step_number[]
     initial_extrap = irkgl_cache.initial_extrap
     len = irkgl_cache.length_u
     tf=irkgl_cache.tf

     s = length(b)
     maxiters = (step_number==1 ? 10+irkgl_cache.maxiters : irkgl_cache.maxiters )
     tj = ttj[1]
     te = ttj[2]
     indices=eachindex(uj)

     dt=dts[1]
     dtprev=dts[2]
     signdt=dts[3]
     sdt=dt*signdt


     step_retcode=true

     if initial_extrap && (step_number>1)
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
    
    
    j_iter = 0  # Counter of fixed_point iterations
    nf=0
    Dmin .= Inf

    iter = true # Initialize iter outside the for loop
    plusIt= true
    diffU = true

    
    @inbounds while (j_iter<maxiters && iter)  
    
        iter = false
          j_iter += 1
    
          for is in 1:s
              nf+=1
              f(F[is], U[is], p,  muladd(sdt,c[is], tj))
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
          
    end # while
  
    if iter  # iter=false implies that j_iter==maxiters
        
        @warn "Interrupted. Reached maximum number of iterations (maxiters=$maxiters). The value dt=$dt may be too large."
        step_retcode=false
     
    end
    
    if  step_retcode

        @inbounds if (j_iter<maxiters && diffU)   
        
              j_iter += 1
    
              for is in 1:s
                  nf+=1
                  f(F[is], U[is], p,  muladd(sdt,c[is],tj))
                  for k in indices
                      L[is][k] = sdt*(b[is]*F[is][k])
                  end
              end
        end
    
    
        @inbounds for k in indices    #Equivalent to compensated summation
    
          L_sum = L[1][k]
          for is in 2:s
              L_sum+=L[is][k]
          end
          res = Base.TwicePrecision(uj[k], ej[k]) + L_sum
    
          uj[k] = res.hi
          ej[k] = res.lo
        end
    
        dtmax = abs((tf-ttj[1])-ttj[2])
        if abs(sdt) >= dtmax
            ttj[1] = tf
            ttj[2] = 0
        else
            res = Base.TwicePrecision(tj, te) + sdt
            ttj[1] = res.hi
            ttj[2] = res.lo
        end
    
        dts[1]=min(abs(sdt),dtmax)
        dts[2]=dt

        stats.nfpiter+=j_iter
        stats.nf+=nf

    end

    return step_retcode
    
end
    

function IRKGLstep_adap!(ttj,uj,ej,dts,stats,irkgl_cache::NbodyIRKGL16.IRKGL_Cache{uT,tT,fT,pT}) where {uT,tT,fT,pT}

    f = irkgl_cache.odef
    p = irkgl_cache.p
    b = irkgl_cache.b
    c = irkgl_cache.c
    d = irkgl_cache.d
    a = irkgl_cache.a
    mu = irkgl_cache.mu
    nu = irkgl_cache.nu
    theta=irkgl_cache.theta
    omega=irkgl_cache.omega
    alpha=irkgl_cache.alpha
    K = irkgl_cache.K
    logK = irkgl_cache.logK
    Kinv = irkgl_cache.Kinv
    Tau = irkgl_cache.Tau
    Tau_ = irkgl_cache.Tau_
    U = irkgl_cache.U
    U_ = irkgl_cache.U_
    L = irkgl_cache.L
    F = irkgl_cache.F
    Dmin = irkgl_cache.Dmin
    step_number = irkgl_cache.step_number[]
    initial_extrap = irkgl_cache.initial_extrap
    len = irkgl_cache.length_u
    tf=irkgl_cache.tf
    Dtau=irkgl_cache.Dtau

    R=tT(2)^10 
    
    s = length(b)
    extra_iters = (step_number>1 && initial_extrap ? 1 : 9)
    maxiters = irkgl_cache.maxiters + extra_iters - 1
    tj = ttj[1]
    te = ttj[2] 
    dtmax = abs((tf-tj)-te) 
    indices=eachindex(uj)
   
    dtaux = dts[1]
    dtaux_ = dtaux
    dt = min(dtaux,dtmax)      # Initialized as dt=0 at first step
    dtprev=dts[2]
    signdt=dts[3]
    sdt=signdt*dt
    
    step_retcode = true
    nf=0
    
    if initial_extrap && step_number>1
        
          for is in 1:s
             for k in indices
                 dUik = muladd(nu[is,1], L[1][k], ej[k])
                 for js in 2:s
                    dUik = muladd(nu[is,js], L[js][k], dUik)
                 end
                 U_[is][k] =  dUik  # U = uj + U_
             end
          end
        
          lambda=dt/dtprev-1
       
          if abs(lambda)>eps(tT) 
            
               for i in 1:s
                   sumbetai=zero(tT)
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
                            dUik = muladd(alpha[is,js+1], U_[js][k], dUik)
                        end
                        U[is][k] = uj[k] + dUik
                   end
               end              
            
          else    
               for is in 1:s
                   for k in indices
                        U[is][k] = uj[k] + U_[is][k]
                   end
               end 
            
          end  

          aux=zero(tT)
          daux=zero(tT)
          for is in 1:s
            nf+=1
            kf=f(F[is], U[is], p,  muladd(sdt,c[is],tj))
            aux = muladd(b[is],kf,aux)
            daux = muladd(d[is],kf,daux)
          end       
                 
          diff_dt = (Dtau - dt*aux)/daux
          dtaux = NbodyIRKGL16.Rdigits(dt + diff_dt,R)
        
    else
          for is in 1:s
              for k in indices
                   U[is][k] = uj[k] + ej[k]
              end
          end
          aux=zero(tT)
          for is in 1:s
            nf+=1
            kf=f(F[is], U[is], p,  muladd(sdt,c[is],tj))
            aux = muladd(b[is],kf,aux)
          end
          
          dtaux = NbodyIRKGL16.Rdigits(Dtau/aux,R)

    end
  

    iter = true # Initialize iter outside the for loop
    plusIt= true
    diffU = true
    j_iter = 1
    
    if ((abs(dtaux-dtaux_)>0.5*max(dtaux,dtaux_)) && (step_number>1)) || (dtaux<=0) 

        if dtaux<=0
             @warn "dtaux<=0 --> dtaux=dtaux_, n=$step_number, j=$j_iter,  dt_=$dtaux_, dt=$dtaux"
        end
        
        @debug("n=$step_number, j=$j_iter, tj=$(Float32(tj)), dtaux_=$(Float32(dtaux_)), dtaux=$(Float32(dtaux)), dt_re=$(Float64((abs(dtaux-dtaux_)/max(dtaux,dtaux_))))")
        if initial_extrap && (step_number>1)
            for is in 1:s
                for k in indices
                    U[is][k] = uj[k] + ej[k]
                end
            end

            aux=zero(tT)
            for is in 1:s
                nf+=1
                kf=f(F[is], U[is], p,  muladd(sdt,c[is],tj))
                aux = muladd(b[is],kf,aux)
            end

            initial_extrap = false
            dtaux=NbodyIRKGL16.Rdigits(Dtau/aux,R)
            Dmin .= Inf
            diffU = true
            j_iter = 1
        else
            @warn "Convergence failure. Numerical integration interrupted. \n The value of dt=$Dtau appears to be too large. Please try again with a smaller value of dt"
            iter = false 
            step_retcode=false
        end
    end      

    dt_new=min(dtaux,dtmax)
    sdt_new=signdt*dt_new  
    dsdt = sdt_new-sdt

    Dmin .= Inf
   
    @inbounds while (j_iter<maxiters && iter)  
        
        @debug("n=$step_number, j=$j_iter, tj=$(Float32(tj)), dtaux_=$(Float32(dtaux_)), dtaux=$(Float32(dtaux)), dt_re=$(Float64((abs(dtaux-dtaux_)/max(dtaux,dtaux_))))")
        iter = false
         
        for is in 1:s
            for k in indices
                U_[is][k] = U[is][k]
            end
        end
                     
        for is in 1:s
             for js in 1:s
                 mu[js,is] = sdt*a[js,is] 
              end
              mu[is,is] += c[is]*dsdt
        end
        
        for is in 1:s
             for k in indices
                dUik = ej[k]
                for js in 1:s
                   dUik = muladd(mu[is,js], F[js][k], dUik)
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


        if iter

            dt = dt_new
            sdt = sdt_new

            aux=zero(tT)
            daux=zero(tT)
            for is in 1:s
                nf+=1
                kf=f(F[is], U[is], p,  muladd(sdt,c[is],tj))
                aux = muladd(b[is],kf,aux)
                daux = muladd(d[is],kf,daux)
                K[is] = kf
            end

            j_iter += 1
            diff_dt = (Dtau - dt*aux)/daux 

            if abs(diff_dt)<=eps(dt) #|| true
                diff_dt = zero(tT)
            end

            dtaux_ = dtaux
            dtaux = NbodyIRKGL16.Rdigits(dt + diff_dt,R)  

            if (abs(dtaux-dtaux_)>0.5*max(dtaux,dtaux_)) || (dtaux<=0)

                if dtaux<=0
                     @warn "(while) dtaux<=0 --> dtaux=Dtau/aux, n=$step_number, j=$j_iter, dt_=$dtaux_, dt=$dtaux"
                end
                
                @debug("n=$step_number, j=$j_iter, tj=$(Float32(tj)), dtaux_=$(Float32(dtaux_)), dtaux=$(Float32(dtaux)), dt_re=$(Float64((abs(dtaux-dtaux_)/max(dtaux,dtaux_))))")
                
                if initial_extrap && (step_number>1)
                    for is in 1:s
                        for k in indices
                            U[is][k] = uj[k] + ej[k]
                        end
                    end

                    aux=zero(tT)
                    for is in 1:s
                        nf+=1
                        kf=f(F[is], U[is], p,  muladd(sdt,c[is],tj))
                        aux = muladd(b[is],kf,aux)

                    end

                    initial_extrap = false
                    dtaux = NbodyIRKGL16.Rdigits(Dtau/aux,R)

                    Dmin .= Inf
                    diffU = true
                    j_iter = 1
                else
                    @warn "Convergence failure. Numerical integration interrupted. \n The value of dt=$Dtau appears to be too large. Please try again with a smaller value of dt"
                    iter = false
                    step_retcode=false
                end
            end      
            
            dt_new = min(dtaux, dtmax )
            sdt_new = signdt*dt_new  
            dsdt = sdt_new - sdt

        end

    end # while
 
    dt = dt_new
    sdt = sdt_new     
            
    if (iter && (j_iter==maxiters) && step_retcode)
            @warn "Numerical integration interrupted. Reached maximum number of iterations (maxiters=$maxiters). \n The value of dt=$Dtau appears to be too large. Please try again with a smaller value of dt"
            step_retcode=false
    end   
    

    if  step_retcode
    
        @inbounds if (j_iter<maxiters && diffU)   
        
                j_iter += 1
        
                for is in 1:s

                    nf += 1
                    kf = f(F[is], U[is], p,  muladd(sdt,c[is],tj))
                    K[is] = kf
                  
                end

        end 

        
        for is in 1:s    
            for k in indices
                L[is][k] = sdt*(b[is]*F[is][k])
            end
        end


        @inbounds for k in indices    #Equivalent to compensated summation
    
            L_sum = L[1][k]
            for is in 2:s
                L_sum += L[is][k]
            end
            res = Base.TwicePrecision(uj[k], ej[k]) + L_sum
    
            uj[k] = res.hi
            ej[k] = res.lo
        end
    
        if dtaux >= dtmax
                ttj[1] = tf
                ttj[2] = 0
        else
                res = Base.TwicePrecision(tj, te) + sdt
                ttj[1] = res.hi
                ttj[2] = res.lo
        end

        @. mu = dt * a     

        for is in 1:s
                    Tau_i = zero(tT)
                    for js in 1:s
                    Tau_i = muladd(mu[is,js], K[js], Tau_i)
                    end
                    Tau[is] =  Tau_i
                    logK[is] = log(K[is])
        end 

        
        NbodyIRKGL16.PolInterp!(K, Kinv, Tau, logK, s, Tau_)  # Tau_ = (1+c)*Dtau
        @. Kinv = exp(-K)
            
        dtaux = zero(tT)
        for js in 1:s
            dtaux =muladd(b[js], Kinv[js], dtaux)
        end
        dtaux = NbodyIRKGL16.Rdigits(Dtau*dtaux,R)

        dts[1]=dtaux  # A guess for next time-step
        dts[2]=dt
            
        stats.nfpiter+=j_iter
        stats.nf+=nf


    end
        
    return step_retcode 

   
end
