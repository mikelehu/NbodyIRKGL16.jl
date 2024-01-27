
#
#  launch_IRKGL16_tests: SCIML IRKGL16-bertsioa
#  launch_IRKGL_tests: 2022. urtean garatutako IRKGL
#  launch_method_tests:  SCIML integratzaileak


struct WPTests
   errors::Array{Float64,1}
   times::Array{Float64,1}
end



function launch_IRKGL16_tests(method, final_state, prob, launch_list; adaptive=true, initial_interp=true, nruns=10)

#
#    @belapsed erabiltzea problemak ematen ditu!!!
#
     k=length(launch_list)
     errors=zeros(k)
     times=zeros(k)

     @unpack f,u0,tspan,p,kwargs=prob
     t0=tspan[1]
     tf=tspan[2]

     if (adaptive==true)

      tols=launch_list

      for i in 1:k

        tol=tols[i]

        soli=solve(prob,method, abstol=tol, reltol=tol, adaptive=true, save_everystep=false)
        errors[i]=norm(final_state-soli.u[end])/norm(final_state)

        for k in 1:nruns
           times[i]+=@elapsed solve(prob,method,  abstol=tol, reltol=tol, adaptive=true, save_everystep=false)
        end

        times[i]=times[i]/nruns

      end
    
    else #adaptive==false

      dts=launch_list

      for i in 1:k

        dti=dts[i]
         #n=Int64((tf-t0)/dti)
         # m=n => save_everystep=false

       soli=solve(prob,method, dt=dti, adaptive=false, save_everystep=false)
       errors[i]=norm(final_state-soli.u[end])/norm(final_state)

       for k in 1:nruns
          times[i]+=@elapsed solve(prob,method, dt=dti, adaptive=false, save_everystep=false)
       end

       times[i]=times[i]/nruns

     end
    
    end


     WPTests(errors,times)

end




function launch_IRKGL_tests(integrator, final_state, prob, s, dts; dim=1, adaptive=false, initial_interp=1, floatType=Float64, maxiters=100, nruns=10)

  #
  #    @belapsed erabiltzea problemak ematen ditu!!!
  #
  #    integrator= IRKGL_seq, IRKGL_simd
  #                IRKNGL_seq, IRKNGL_simd
  #                IRKNGL_seq, IRKNGL_simd
  #

       k=length(dts)
       errors=zeros(k)
       times=zeros(k)  
       
       if (integrator==IRKGL_simd)   alg=IRKGL_simd(s=s, initial_interp=initial_interp, dim=dim, floatType=floatType) end
       if (integrator==IRKNGL_simd)  alg=IRKNGL_simd(s=s, initial_interp=initial_interp, dim=dim, floatType=floatType) end
       if (integrator==IRKGL_Seq)    alg=IRKGL_Seq(s=s, initial_interp=initial_interp) end
       if (integrator==IRKNGL_Seq)   alg=IRKNGL_Seq(s=s, initial_interp=initial_interp) end
      

       if (adaptive==true)

        for i in 1:k
  
          dti=dts[i]
          soli=solve(prob,alg; dt=dti, adaptive=true,  save_everystep=false, maxiters=maxiters)
          errors[i]=norm(final_state-soli.u[end])/norm(final_state)
  
          for k in 1:nruns
             times[i]+=@elapsed solve(prob,alg; dt=dti, adaptive=true, save_everystep=false, maxiters=maxiters)
          end
  
          times[i]=times[i]/nruns
  
        end

       
       else  # adaptive==false

        for i in 1:k
  
          dti=dts[i]
          soli=solve(prob,alg; dt=dti, adaptive=false,  save_everystep=false, maxiters=maxiters)
          errors[i]=norm(final_state-soli.u[end])/norm(final_state)
  
          for k in 1:nruns
             times[i]+=@elapsed solve(prob,alg; dt=dti, adaptive=false, save_everystep=false, maxiters=maxiters)
          end
  
          times[i]=times[i]/nruns
  
        end

      end
  
       WPTests(errors,times)
  
end



function launch_method_tests(method, final_state, prob, launch_list; adaptive=true, maxiters=10^9,  nruns=10)

#
#    @belapsed erabiltzea problemak ematen ditu!!!
#

     k=length(launch_list)
     errors=zeros(k)
     times=zeros(k)

     if (adaptive==true)

       tols=launch_list

       for i in 1:k

           tol=tols[i]
           soli= solve(prob, method, abstol=tol, reltol=tol, adaptive=true, save_everystep=false, dense=false, maxiters=maxiters);
           errors[i]=norm(final_state-soli.u[end])/norm(final_state)

           for k in 1:nruns
               times[i]+=@elapsed solve(prob, method, abstol=tol, reltol=tol, adaptive=true, save_everystep=false, dense=false, maxiters=maxiters);
           end

           times[i]=times[i]/nruns

     end

    else # adaptive_false

      dts=launch_list

      for i in 1:k

          dti=dts[i]
          soli= solve(prob, method, dt=dti, adaptive=false, save_everystep=false, dense=false)
          errors[i]=norm(final_state-soli.u[end])/norm(final_state)

          for k in 1:nruns
              times[i]+=@elapsed solve(prob, method, dt=dti, adaptive=false, save_everystep=false, dense=false);
          end

          times[i]=times[i]/nruns

      end

    end

    WPTests(errors,times)

end
