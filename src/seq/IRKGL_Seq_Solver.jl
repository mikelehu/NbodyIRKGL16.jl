@enum retcode begin
    Success
    Failure
end

struct IRKGL_Cache{uT,tT,fT,pT}
    odef::fT # function defining the ODE system
    p::pT # parameters and so
    b::Vector{tT}
    c::Vector{tT}
    mu::Array{tT,2}
    nu::Array{tT,2}
    d::Vector{tT}
    U::Vector{uT}
    U_::Vector{uT}
    L::Vector{uT}
    F::Vector{uT}
    Dmin::Array{tT,1}
    maxiters::Int64
    step_number::Array{Int64,0}
    initial_extrap::Array{Int64,0}
    length_u::Int64
    tf::tT
    Dtau::tT
    step_retcode::retcode
end



abstract type IRKAlgorithm{s,initial_extrap,m, Dtau} <: OrdinaryDiffEqAlgorithm end
struct IRKGL_seq{s, initial_extrap, m, Dtau} <: IRKAlgorithm{s, initial_extrap, m, Dtau} end
#IRKGL_seq(;s=8, initial_extrap=1, m=1, Dtau=one(eltype(tspanType)))=IRKGL_seq{s, initial_extrap, m, Dtau}()
IRKGL_seq(;s=8, initial_extrap=1, m=1,Dtau=1.)=IRKGL_seq{s, initial_extrap, m, Dtau}()

function DiffEqBase.__solve(prob::DiffEqBase.AbstractODEProblem{uType,tspanType,isinplace},
    alg::IRKGL_seq{s,initial_extrap, m, Dtau}, args...;
#    dt=0.,
    dt=zero(eltype(tspanType)),
    save_everystep=true,
    adaptive=false,
    maxiters=100,
    kwargs...) where {uType,tspanType,isinplace,s, initial_extrap, m, Dtau}



    trace = false
    checks=true

    stats = DiffEqBase.Stats(0)
    stats.nf=0
    stats.nnonliniter=0
    stats.naccept=0
    
    @unpack f,u0,tspan,p,kwargs=prob
    f= SciMLBase.unwrapped_f(prob.f) 


    uiType=eltype(uType)
    if uiType <:Complex 
        uiType=real(uiType)
    end

    tType=eltype(tspanType)
    Dtau0=convert(tType,Dtau)

#   Memory preallocation (IRKL_Cache)   

    step_number = Array{Int64,0}(undef)
    step_number[] = 0
    init_interp =  Array{Int64,0}(undef)
    init_interp[] = initial_extrap

    (b, c, mu, nu, d) = IRKGLCoefficients(s,dt)
    length_u = length(u0)
    indices=1:length_u
    step_retcode::retcode=Success


    U=Array{uType}(undef, s)
    for i in 1:s
         U[i]=zero(u0)
    end

    U_=deepcopy(U)
    L=deepcopy(U)
    F=deepcopy(U)

    Dmin=Array{uiType}(undef,length_u)
    Dmin.=zero(uiType)
 

    dtprev=zero(tType)
    signdt=sign(tspan[2]-tspan[1]) 
           
    if signdt==1          # forward integration
        t0=prob.tspan[1]
        tf=prob.tspan[2]   
    else                  # backward integration
        t0=prob.tspan[2]
        tf=prob.tspan[1]   
    end

    dts=[dt,dtprev,signdt]

    irkgl_cache = IRKGL_Cache(f,p,b,c,mu,nu,d,
                              U, U_, L, 
                              F,Dmin,maxiters,step_number,
                              init_interp,length_u,tf, Dtau0, step_retcode)
#

    uu = uType[]
    tt = tType[]
                                
    push!(uu,copy(u0))
    push!(tt,t0)
                                                           
    tj = [t0, zero(t0)]
    uj = copy(u0)
    ej= zero(u0)

    if adaptive==false && dt==0.
        println("Error: dt is required.")
        checks=false
    end


#    if adaptive==true
#        println("Error: only allowed adaptive =false")
#        checks=false
#    end
     
    if checks==true

        cont=true
        error_warn=0

        if adaptive
        
            while cont

                for i in 1:m 

                    step_number[]+= 1
                    IRKGLstep_adap!(tj,uj,ej,dts,stats,irkgl_cache)
               
                    if (step_retcode==Failure)
                        error_warn=1
                        break
                    end
                    
                    if (dts[1]==0)  cont=false  end

                end

                if save_everystep !=false
                    push!(uu,uj+ej)
                    push!(tt,tj[1]+tj[2])
                end

            end # end while
        
        else

            while cont

                for i in 1:m 

                    step_number[]+= 1
                    IRKGLstep_fixed!(tj,uj,ej,dts,stats,irkgl_cache)
               
                    if (step_retcode==Failure)
                        error_warn=1
                        break
                    end
                    
                    if (dts[1]==0)  cont=false  end

                end

                if save_everystep !=false
                    push!(uu,uj+ej)
                    push!(tt,tj[1]+tj[2])
                end

            end # end while

        end

        stats.naccept=step_number[]
        if error_warn!=0

            println("Error during the integration warn=$error_warn")
            sol=DiffEqBase.build_solution(prob,alg,tt,uu,stats=stats,retcode= ReturnCode.Failure)
        
        else
            if tt[end]!=tf
                push!(uu,uj+ej)
                push!(tt,tj[1]+tj[2])
            end

            sol=DiffEqBase.build_solution(prob,alg,tt,uu,stats=stats,retcode= ReturnCode.Success)
        end

        return(sol)


    end
  
end


