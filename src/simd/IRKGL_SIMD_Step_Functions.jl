#   
#  IRKGL Step functions

#      IRKGLstep_SIMD_fixed!   
#      IRKGLstep_SIMD_adap!   


function IRKGLstep_SIMD_fixed!(ttj, uj,ej,dts,stats,irkglcache::IRKGL_SIMD.IRKGL_SIMD_Cache{floatT,fType,pType,s_,dim_}) where {floatT,fType,pType,s_,dim_}


    f = irkglcache.odef
    p = irkglcache.p
    b = irkglcache.b
    c = irkglcache.c
    mu = irkglcache.mu
    nu = irkglcache.nu
    U = irkglcache.U
    U_ = irkglcache.U_
    L = irkglcache.L
    F = irkglcache.F
    Dmin = irkglcache.Dmin
    step_number = irkglcache.step_number[]
    initial_extrap = irkglcache.initial_extrap
    len = irkglcache.length_u
    tf= irkglcache.tf

    s = length(b)
    maxiters = (step_number==1 ? 10+irkglcache.maxiters : irkglcache.maxiters )
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
         Lk = IRKGL_SIMD.getindex_(L,k)
         dUk = muladd(nu[1], Lk[1], ej[k])
         for js in 2:s
             dUk = muladd(nu[js], Lk[js], dUk)
         end
         IRKGL_SIMD.setindex_!(U, uj[k]+dUk, k)
      end

    else

       for k in indices
         uej = uj[k] + ej[k]
         IRKGL_SIMD.setindex_!(U, uej, k)
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

             Fk = IRKGL_SIMD.getindex_(F,k)
             Lk = sdt*(b*Fk)
             dUk = muladd(mu[1], Lk[1], ej[k])
             for is in 2:s
                     dUk = muladd(mu[is], Lk[is], dUk)
             end
             Uk = uj[k]+dUk
             IRKGL_SIMD.setindex_!(U, Uk, k)
             IRKGL_SIMD.setindex_!(L, Lk, k)
             Uk_ = IRKGL_SIMD.getindex_(U_,k)
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
                    Fk = IRKGL_SIMD.getindex_(F, k)
                    Lk = sdt*(b*Fk)
                    IRKGL_SIMD.setindex_!(L, Lk, k)
                end
        end


        @inbounds for k in indices    #Equivalent to compensated summation
            Lk = IRKGL_SIMD.getindex_(L,k)
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



function IRKGLstep_SIMD_adap!(ttj, uj,ej,dts,stats,irkglcache::IRKGL_SIMD.IRKGL_SIMD_Cache{tT,fType,pType,s_,dim_}) where {tT,fType,pType,s_,dim_}


    f = irkglcache.odef
    p = irkglcache.p
    b = irkglcache.b
    c = irkglcache.c
    d = irkglcache.d
    a = irkglcache.a
    s_beta=irkglcache.s_beta   # beharrezkoa da 2024-04-24
    mu = irkglcache.mu
    nu = irkglcache.nu
    theta = irkglcache.theta
    omega = irkglcache.omega
    alpha = irkglcache.alpha
    K = irkglcache.K
    logK = irkglcache.logK
    Kinv = irkglcache.Kinv
    Tau = irkglcache.Tau
    Tau_ = irkglcache.Tau_
    U = irkglcache.U
    U_ = irkglcache.U_
    L = irkglcache.L
    F = irkglcache.F
    Dmin = irkglcache.Dmin
    step_number = irkglcache.step_number[]
    initial_extrap = irkglcache.initial_extrap
    len = irkglcache.length_u
    tf= irkglcache.tf
    Dtau = irkglcache.Dtau

    R=tT(2)^10 

    s = length(b)
    extra_iters = (step_number>1 && initial_extrap ? 1 : 9)
    maxiters = irkglcache.maxiters + extra_iters - 1
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
            Lk = IRKGL_SIMD.getindex_(L,k)
            dUk = muladd(nu[1], Lk[1], ej[k])
            for js in 2:s
                dUk = muladd(nu[js], Lk[js], dUk)
            end
            IRKGL_SIMD.setindex_!(U_, dUk, k)       # U = uj + U_
        end

        lambda=dt/dtprev-1

        if abs(lambda)>eps(tT) 

            s_beta=0
            
            for j in 1:s+1

                thetaj=IRKGL_SIMD.getindex_(theta,j)
                omegaj=IRKGL_SIMD.getindex_(omega,j)
                betaj=1/muladd(lambda,thetaj,omegaj)
                IRKGL_SIMD.setindex_!(alpha, betaj,j)
                s_beta+=betaj

            end

            for j in 1:s+1
               alphaj=IRKGL_SIMD.getindex_(alpha,j)
               IRKGL_SIMD.setindex_!(alpha, alphaj/s_beta, j)
            end              
       
            for k in indices
                dUk=IRKGL_SIMD.getindex_(U_,k)
                dUik = alpha[1]*ej[k]
                for js in 1:s
                    dUik = muladd(alpha[js+1], dUk[js], dUik)
                end
                IRKGL_SIMD.setindex_!(U, uj[k]+dUik, k)
            end


        else    

            for k in indices
                dUk=IRKGL_SIMD.getindex_(U_,k)
                IRKGL_SIMD.setindex_!(U, uj[k]+dUk, k)
            end
     
        end    
   
        nf+=s
        kf=f(F, U, p, muladd(sdt,c,tj))
        aux=sum(b*kf)
        daux=sum(d*kf)
   
        diff_dt = (Dtau - dt*aux)/daux
        dtaux = IRKGL_SIMD.Rdigits(dt + diff_dt,R)

    else

       for k in indices
         uej = uj[k] + ej[k]
         IRKGL_SIMD.setindex_!(U, uej, k)
       end

       nf+=s
       kf=f(F, U, p, muladd(sdt,c,tj))
       aux=sum(b*kf)
       
       dtaux = IRKGL_SIMD.Rdigits(Dtau/aux,R)

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
               IRKGL_SIMD.setindex_!(U, uj[k] + ej[k], k)
            end

            nf+=s
            kf=f(F, U, p, muladd(sdt,c,tj))
            aux=sum(b*kf)

            initial_extrap = false
            dtaux = IRKGL_SIMD.Rdigits(Dtau/aux,R)
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
            aj=IRKGL_SIMD.getindex_(a,j)
            IRKGL_SIMD.setindex_!(mu,sdt*aj,j)
            muj=IRKGL_SIMD.getindex_(mu,j)
            mujj=muj[j]+c[j]*dsdt
            mu.data[j,j]=mujj
        end
        

        diffU = false

        for k in indices

             Fk = IRKGL_SIMD.getindex_(F,k)
             dUk = muladd(mu[1], Fk[1], ej[k])
             for is in 2:s
                     dUk = muladd(mu[is], Fk[is], dUk)
             end
             Uk = uj[k]+dUk

             IRKGL_SIMD.setindex_!(U, Uk, k)
             Uk_ = IRKGL_SIMD.getindex_(U_,k)
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
            dtaux = IRKGL_SIMD.Rdigits(dt + diff_dt,R)  
           
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
                    dtaux = IRKGL_SIMD.Rdigits(Dtau/aux,R)
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

            Fk = IRKGL_SIMD.getindex_(F, k)
            Lk = sdt*(b*Fk)
            IRKGL_SIMD.setindex_!(L, Lk, k)
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
            aj=IRKGL_SIMD.getindex_(a,j)
            IRKGL_SIMD.setindex_!(mu,dt*aj,j)
        end

        Tau=muladd(mu[1],K[1],0)
        for j in 2:s
            Tau=muladd(mu[j],K[j],Tau)
        end

        logK=log(K)
    
        K=IRKGL_SIMD.PolInterp!(K, Kinv, Tau, logK, s, Tau_)  # Tau_ = (1+c)*Dtau
        #for i in 1:s
        #    pz=IRKGL_SIMD.PolInterp(Tau,logK,s,Tau_[i])
        #    K=Base.setindex(K,pz,i)
        #end

        Kinv = exp(-K)

        #if sum(Kinv)==Inf

        #    @warn "Numerical integration interrupted. \n The value of dt=$Dtau appears to be too small. Please try again with a greater value of dt"
        #    step_retcode=false

        #else

        dtaux= IRKGL_SIMD.Rdigits(Dtau*sum(b*Kinv),R)

        dts[1]=dtaux  # A guess for next time-step
        dts[2]=dt

        stats.nnonliniter+=j_iter
        stats.nf+=nf
        
        #end

    end

    return step_retcode

end
