#
#  IRKGLstep_fixed!
#  IRKGLstep_adap!


function IRKGLstep_fixed!(ttj,uj,ej,dts,stats,irkgl_cache::IRKGL_Cache{uT,tT,fT,pT}) where {uT,tT,fT,pT}

    
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
     initial_extrap = irkgl_cache.initial_extrap[]
     len = irkgl_cache.length_u
     tf=irkgl_cache.tf
     step_retcode = irkgl_cache.step_retcode


     s = length(b)
     maxiters = (step_number==1 ? 10+irkgl_cache.maxiters : irkgl_cache.maxiters )
     tj = ttj[1]
     te = ttj[2]
     indices=1:len

     dt=dts[1]
     dtprev=dts[2]
     signdt=dts[3]
     sdt=dt*signdt
    

     if (initial_extrap==1)
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

    
      end # while
    
    if (!iter && (j_iter==maxiters))
         println("Failure !!!  Step=",tj+te, " dt=", dts[1])
         step_retcode=Failure
         return nothing
      
    else 
    
        @inbounds if (j_iter<maxiters && diffU)   
        
              j_iter += 1
    
              for is in 1:s
                  nf+=1
                  f(F[is], U[is], p,  tj  + sdt*c[is])
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
    
    
        res = Base.TwicePrecision(tj, te) + sdt
        ttj[1] = res.hi
        ttj[2] = res.lo

    
        dts[1]=min(abs(dt),abs(tf-(ttj[1]+ttj[2])))
        dts[2]=dt

        if dts[2]!=Inf
            stats.nnonliniter+=j_iter
            stats.nf+=nf
        end

        retcode=Success
        return nothing

    end

    
end
    
function IRKGLstep_adap!(ttj,uj,ej,dts,stats,irkgl_cache::IRKGL_Cache{uT,tT,fT,pT}) where {uT,tT,fT,pT}

    
    f = irkgl_cache.odef
    p = irkgl_cache.p
    b = irkgl_cache.b
    c = irkgl_cache.c
#    d = irkgl_cache.d
    d = irkgl_cache.b
    mu = irkgl_cache.mu
    nu = irkgl_cache.nu
    U = irkgl_cache.U
    U_ = irkgl_cache.U_
    L = irkgl_cache.L
    F = irkgl_cache.F
    Dmin = irkgl_cache.Dmin
    step_number = irkgl_cache.step_number[]
    initial_extrap = irkgl_cache.initial_extrap[]
    len = irkgl_cache.length_u
    tf=irkgl_cache.tf
    Dtau=irkgl_cache.Dtau
    step_retcode = irkgl_cache.step_retcode


    s = length(b)
    maxiters = (step_number==1 ? 10+irkgl_cache.maxiters : irkgl_cache.maxiters )
    tj = ttj[1]
    te = ttj[2]
    indices=1:len
   
    dtprev=dts[2]
    signdt=dts[3]

#    if (initial_extrap==1)
#        for is in 1:s
#            for k in indices
#                dUik = muladd(nu[is,1], L[1][k], ej[k])
#                for js in 2:s
#                   dUik = muladd(nu[is,js], L[js][k], dUik)
#                end
#                U[is][k] =  uj[k]  + dUik
#            end
#        end
#    else
        for is in 1:s
            for k in indices
               U[is][k] = uj[k] + ej[k]
            end
        end
#    end

    nf=0
    sdt=zero(eltype(tj))
    aux=zero(eltype(uj))
    for is in 1:s
        nf+=1
        kf=f(F[is], U[is], p,  tj  + sdt*c[is])
        aux+=d[is]*kf
    end


    dt=min(abs(Dtau/aux),abs(tf-(ttj[1]+ttj[2]) ))
    sdt=signdt*dt  
   
   j_iter = 0  # Counter of fixed_point iterations
   Dmin .= Inf

   iter = true # Initialize iter outside the for loop
   plusIt= true
   diffU = true

#   println("urratsa:", tj)
   
   @inbounds while (j_iter<maxiters && iter)  
   
         iter = false
         j_iter += 1

 #        println("iterazioa. j=",j_iter, ",sdt=", sdt )
         

        for is in 1:s
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

        if iter
            aux=zero(eltype(uj))
            for is in 1:s
                nf+=1
                kf=f(F[is], U[is], p,  tj  + sdt*c[is])
                aux+=d[is]*kf
            end
            dt=min(abs(Dtau/aux),abs(tf-(ttj[1]+ttj[2]) ))
            sdt=signdt*dt
        end
   
     end # while
   
   if (!iter && (j_iter==maxiters))
        println("Failure !!!  Step=",tj+te, " dt=", dts[1])
        step_retcode=Failure
        return nothing
     
   else 
   
       @inbounds if (j_iter<maxiters && diffU)   
       
             j_iter += 1
   
             for is in 1:s
                 nf+=1
                 f(F[is], U[is], p,  tj  + sdt*c[is])
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
   
   
       res = Base.TwicePrecision(tj, te) + sdt
       ttj[1] = res.hi
       ttj[2] = res.lo

   
       dts[1]=min(abs(dt),abs(tf-(ttj[1]+ttj[2])))
       dts[2]=dt

       if dts[2]!=Inf
           stats.nnonliniter+=j_iter
           stats.nf+=nf
       end

       retcode=Success
       return nothing

   end

   
end