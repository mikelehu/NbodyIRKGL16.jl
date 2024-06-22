#   
#  IRKGL Step functions

#      IRKGLstep_SIMD_fixed!   
#      IRKGLstep_SIMD_adap!   


function IRKGLstep_SIMD_fixed!(ttj, uj,ej,dts,stats,irkgl_cache::NbodyIRKGL16.IRKGL_SIMD_Cache{floatT,fType,pType,s_,dim_}) where {floatT,fType,pType,s_,dim_}


    f = irkgl_cache.odef
    p = irkgl_cache.p
    b = irkgl_cache.b
    c = irkgl_cache.c
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
    tf= irkgl_cache.tf

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

      for k in indices
         Lk = NbodyIRKGL16.getindex_(L,k)
         dUk = muladd(nu[1], Lk[1], ej[k])
         for js in 2:s
             dUk = muladd(nu[js], Lk[js], dUk)
         end
         NbodyIRKGL16.setindex_!(U, uj[k]+dUk, k)
      end

    else

       for k in indices
         uej = uj[k] + ej[k]
         NbodyIRKGL16.setindex_!(U, uej, k)
       end

    end


    j_iter = 0  # counter of fixed_point iterations
    nf=0
    Dmin .= Inf

    iter = true # Initialize iter outside the for loop
    plusIt=true
    diffU = true

    @inbounds while (j_iter<maxiters && iter ) 

         iter = false
         j_iter += 1

         U_.data .= U.data
         nf+=s
         f(F, U, p, muladd(sdt,c,tj))

         diffU = false

         for k in indices

             Fk = NbodyIRKGL16.getindex_(F,k)
             Lk = sdt*(b*Fk)
             dUk = muladd(mu[1], Lk[1], ej[k])
             for is in 2:s
                     dUk = muladd(mu[is], Lk[is], dUk)
             end
             Uk = uj[k]+dUk
             NbodyIRKGL16.setindex_!(U, Uk, k)
             NbodyIRKGL16.setindex_!(L, Lk, k)
             Uk_ = NbodyIRKGL16.getindex_(U_,k)
             DY = maximum(abs(Uk-Uk_))

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

    if  iter # iter=false implies that j_iter==maxiters

        @warn "Interrupted. Reached maximum number of iterations (maxiters=$maxiters). The value dt=$dt may be too large."
        step_retcode=false

    end

    if  step_retcode   

        @inbounds if (j_iter<maxiters && diffU) 
                j_iter += 1
                nf+=s
                f(F, U, p, muladd(sdt,c,tj))

                for k in indices
                    Fk = NbodyIRKGL16.getindex_(F, k)
                    Lk = sdt*(b*Fk)
                    NbodyIRKGL16.setindex_!(L, Lk, k)
                end
        end


        @inbounds for k in indices    #Equivalent to compensated summation
            Lk = NbodyIRKGL16.getindex_(L,k)
            L_sum = sum(Lk)
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

        stats.nnonliniter+=j_iter
        stats.nf+=nf

    end

    return step_retcode


end



function IRKGLstep_SIMD_adap!(ttj, uj,ej,dts,stats,irkgl_cache::NbodyIRKGL16.IRKGL_SIMD_Cache{tT,fType,pType,s_,dim_}) where {tT,fType,pType,s_,dim_}


    f = irkgl_cache.odef
    p = irkgl_cache.p
    b = irkgl_cache.b
    c = irkgl_cache.c
    d = irkgl_cache.d
    a = irkgl_cache.a
    s_beta=irkgl_cache.s_beta   # beharrezkoa da 2024-04-24
    mu = irkgl_cache.mu
    nu = irkgl_cache.nu
    theta = irkgl_cache.theta
    omega = irkgl_cache.omega
    alpha = irkgl_cache.alpha
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
    tf= irkgl_cache.tf
    Dtau = irkgl_cache.Dtau

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

    step_retcode=true      
    nf=0

    if initial_extrap && (step_number>1)

        for k in indices
            Lk = NbodyIRKGL16.getindex_(L,k)
            dUk = muladd(nu[1], Lk[1], ej[k])
            for js in 2:s
                dUk = muladd(nu[js], Lk[js], dUk)
            end
            NbodyIRKGL16.setindex_!(U_, dUk, k)       # U = uj + U_
        end

        lambda=dt/dtprev-1

        if abs(lambda)>eps(tT) 

            s_beta=0
            
            for j in 1:s+1

                thetaj=NbodyIRKGL16.getindex_(theta,j)
                omegaj=NbodyIRKGL16.getindex_(omega,j)
                betaj=1/muladd(lambda,thetaj,omegaj)
                NbodyIRKGL16.setindex_!(alpha, betaj,j)
                s_beta+=betaj

            end

            for j in 1:s+1
               alphaj=NbodyIRKGL16.getindex_(alpha,j)
               NbodyIRKGL16.setindex_!(alpha, alphaj/s_beta, j)
            end              
       
            for k in indices
                dUk=NbodyIRKGL16.getindex_(U_,k)
                dUik = alpha[1]*ej[k]
                for js in 1:s
                    dUik = muladd(alpha[js+1], dUk[js], dUik)
                end
                NbodyIRKGL16.setindex_!(U, uj[k]+dUik, k)
            end


        else    

            for k in indices
                dUk=NbodyIRKGL16.getindex_(U_,k)
                NbodyIRKGL16.setindex_!(U, uj[k]+dUk, k)
            end
     
        end    
   
        nf+=s
        kf=f(F, U, p, muladd(sdt,c,tj))
        aux=sum(b*kf)
        daux=sum(d*kf)
   
        diff_dt = (Dtau - dt*aux)/daux
        dtaux = NbodyIRKGL16.Rdigits(dt + diff_dt,R)

    else

       for k in indices
         uej = uj[k] + ej[k]
         NbodyIRKGL16.setindex_!(U, uej, k)
       end

       nf+=s
       kf=f(F, U, p, muladd(sdt,c,tj))
       aux=sum(b*kf)
       
       dtaux = NbodyIRKGL16.Rdigits(Dtau/aux,R)

    end
  
    iter = true # Initialize iter outside the for loop
    plusIt=true
    diffU = true
    j_iter = 1  

    if ((abs(dtaux-dtaux_)>0.5*max(dtaux,dtaux_)) && (step_number>1)) || (dtaux<=0) 


        if dtaux<=0
             @warn "dtaux<=0 --> dtaux=dtaux_, n=$step_number, j=$j_iter,  dt_=$dtaux_, dt=$dtaux"
        end
        
        @debug("n=$step_number, j=$j_iter, tj=$(Float32(tj)), dtaux_=$(Float32(dtaux_)), dtaux=$(Float32(dtaux)), dt_re=$(Float64((abs(dtaux-dtaux_)/max(dtaux,dtaux_))))")
       
        if initial_extrap && (step_number>1)

            for k in indices
               NbodyIRKGL16.setindex_!(U, uj[k] + ej[k], k)
            end

            nf+=s
            kf=f(F, U, p, muladd(sdt,c,tj))
            aux=sum(b*kf)

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
  
    dt_new=min(dtaux,dtmax)
    sdt_new=signdt*dt_new  
    dsdt = sdt_new-sdt

    Dmin .= Inf


    @inbounds while (j_iter<maxiters && iter ) 

        @debug("n=$step_number, j=$j_iter, tj=$(Float32(tj)), dtaux_=$(Float32(dtaux_)), dtaux=$(Float32(dtaux)), dt_re=$(Float64((abs(dtaux-dtaux_)/max(dtaux,dtaux_))))")

        iter = false
        U_.data .= U.data

        for j in 1:s
            aj=NbodyIRKGL16.getindex_(a,j)
            NbodyIRKGL16.setindex_!(mu,sdt*aj,j)
            muj=NbodyIRKGL16.getindex_(mu,j)
            mujj=muj[j]+c[j]*dsdt
            mu.data[j,j]=mujj
        end
        

        diffU = false

        for k in indices

             Fk = NbodyIRKGL16.getindex_(F,k)
             dUk = muladd(mu[1], Fk[1], ej[k])
             for is in 2:s
                     dUk = muladd(mu[is], Fk[is], dUk)
             end
             Uk = uj[k]+dUk

             NbodyIRKGL16.setindex_!(U, Uk, k)
             Uk_ = NbodyIRKGL16.getindex_(U_,k)
             DY = maximum(abs(Uk-Uk_))

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

            nf+=s
            kf=f(F, U, p, muladd(sdt,c,tj))
            aux=sum(b*kf)
            daux=sum(d*kf)
            K=kf
            
            j_iter += 1
            diff_dt = (Dtau - dt*aux)/daux


            if abs(diff_dt)<=eps(dt)
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

                    for k in indices
                        uej = uj[k] + ej[k]
                        setindex_!(U, uej, k)
                    end

                    nf+=s
                    kf=f(F, U, p, muladd(sdt,c,tj))
                    aux=sum(b*kf)

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

    #if (j_iter==1)
    #    @warn "Numerical integration interrupted. \n The value of Dtau=$Dtau appears to be too small. Please try again with a greater value of dt"
    #    step_retcode=false
    #end

    if  step_retcode

        @inbounds if (j_iter<maxiters && diffU) 
    
            j_iter += 1
    
            nf+=s
            kf=f(F, U, p, muladd(sdt,c,tj))
            K=kf
        
        end

        @inbounds for k in indices    #Equivalent to compensated summation

            Fk = NbodyIRKGL16.getindex_(F, k)
            Lk = sdt*(b*Fk)
            NbodyIRKGL16.setindex_!(L, Lk, k)
            L_sum = sum(Lk)
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


        for j in 1:s
            aj=NbodyIRKGL16.getindex_(a,j)
            NbodyIRKGL16.setindex_!(mu,dt*aj,j)
        end

        Tau=muladd(mu[1],K[1],0)
        for j in 2:s
            Tau=muladd(mu[j],K[j],Tau)
        end

        logK=log(K)
    
        K=NbodyIRKGL16.PolInterp!(K, Kinv, Tau, logK, s, Tau_)  # Tau_ = (1+c)*Dtau
        #for i in 1:s
        #    pz=NbodyIRKGL16.PolInterp(Tau,logK,s,Tau_[i])
        #    K=Base.setindex(K,pz,i)
        #end

        Kinv = exp(-K)

        #if sum(Kinv)==Inf

        #    @warn "Numerical integration interrupted. \n The value of dt=$Dtau appears to be too small. Please try again with a greater value of dt"
        #    step_retcode=false

        #else

        dtaux= NbodyIRKGL16.Rdigits(Dtau*sum(b*Kinv),R)

        dts[1]=dtaux  # A guess for next time-step
        dts[2]=dt

        stats.nnonliniter+=j_iter
        stats.nf+=nf
        
        #end

    end

    return step_retcode

end
