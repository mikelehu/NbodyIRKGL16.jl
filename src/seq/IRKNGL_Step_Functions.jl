#
#  IRKNGLstep_fixed!


function IRKNGLstep_fixed!(ttj, uj,ej, dts, stats, irkngl_cache::IRKNGL_Cache{uT,tT,fT,pT}) where {uT,tT,fT,pT}

 
 f = irkngl_cache.odef
 p = irkngl_cache.p
 b = irkngl_cache.b
 c = irkngl_cache.c
 mu = irkngl_cache.mu
 nu = irkngl_cache.nu
 U = irkngl_cache.U
 U_ = irkngl_cache.U_
 L = irkngl_cache.L
 F = irkngl_cache.F
 Dmin = irkngl_cache.Dmin
 step_number = irkngl_cache.step_number[]
 initial_extrap = irkngl_cache.initial_extrap[]
 len = length(uj)
 lenq = irkngl_cache.length_q
 tf = irkngl_cache.tf
 step_retcode = irkngl_cache.step_retcode


 s = length(b)
 maxiters = (step_number==1 ? 10+irkngl_cache.maxiters : irkngl_cache.maxiters )
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
     for is in 1:s
         for k in indices2
             dUik = muladd(nu[is,1], L[1][k], ej[k])
             for js in 2:s
                dUik = muladd(nu[is,js], L[js][k], dUik)
             end
             U[is][k] =  uj[k]  + dUik
         end
     end
 else
     for is in 1:s
         for k in indices2
            U[is][k] = uj[k] + ej[k]
         end
     end
 end
    
    
 for is in 1:s
     nf+=1
     f(F[is], U[is], p,  tj  + sdt*c[is], 1 )
     for k in indices1
         L[is][k] = sdt*(b[is]*F[is][k])
     end
 end

 for is in 1:s
    for k in indices1
        dUik = muladd(mu[is,1], L[1][k], ej[k])
        for js in 2:s
            dUik = muladd(mu[is,js], L[js][k], dUik)
        end
        U[is][k] =  uj[k] + dUik
    end
 end


 Dmin .= Inf

 iter = true # Initialize iter outside the for loop
 plusIt=true
 diffU = true

  
 @inbounds while (j_iter<maxiters && iter)  

      iter = false
      j_iter += 1

      for is in 1:s
        nf2+=1
        f(F[is], U[is], p,  tj  + sdt*c[is], 2 )
        for k in indices2
            L[is][k] = sdt*(b[is]*F[is][k])
        end
      end

      for is in 1:s
        for k in indices2
           dUik = muladd(mu[is,1], L[1][k], ej[k])
           for js in 2:s
              dUik = muladd(mu[is,js], L[js][k], dUik)
           end
           U_[is][k] = U[is][k]
           U[is][k] =  uj[k] + dUik
        end
      end
        
        
      for is in 1:s
          nf+=1
          f(F[is], U[is], p,  tj  + sdt*c[is], 1 )
          for k in indices1
              L[is][k] = sdt*(b[is]*F[is][k])
          end
      end

      for is in 1:s
          for k in indices1
             dUik = muladd(mu[is,1], L[1][k], ej[k])
             for js in 2:s
                dUik = muladd(mu[is,js], L[js][k], dUik)
             end
           U_[is][k] = U[is][k]
           U[is][k] =  uj[k] + dUik
          end
      end


      diffU = false

      for k in indices1  

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
  
  @inbounds if (j_iter<maxiters && diffU)    #s=8 j_iter>22
      j_iter += 1
        
      for is in 1:s
        nf2+=1
        f(F[is], U[is], p,  tj  + sdt*c[is], 2 )
        for k in indices2
            L[is][k] = sdt*(b[is]*F[is][k])
        end
      end

      for is in 1:s
        for k in indices2
           dUik = muladd(mu[is,1], L[1][k], ej[k])
           for js in 2:s
              dUik = muladd(mu[is,js], L[js][k], dUik)
           end
           U[is][k] =  uj[k] + dUik
        end
      end
        
        
      for is in 1:s
          nf+=1
          f(F[is], U[is], p,  tj  + sdt*c[is], 1 )
          for k in indices1
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
      stats.nf2+=nf2
   end
   
   retcode=Success
   return nothing

  end


end


