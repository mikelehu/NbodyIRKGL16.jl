
struct IRKGL_Cache{uT,tT,fT,pT}
    odef::fT # function defining the ODE system
    p::pT # parameters and so
    b::Vector{tT}
    c::Vector{tT}
    a::Array{tT,2}
    mu::Array{tT,2}
    nu::Array{tT,2}
    theta::Array{tT,2}
    omega::Array{tT,2}
    alpha::Array{tT,2}
    d::Vector{tT}
    K::Vector{tT}
    logK::Vector{tT}
    Kinv::Vector{tT}
    Tau::Vector{tT}
    Tau_::Vector{tT}
    U::Vector{uT}
    U_::Vector{uT}
    L::Vector{uT}
    F::Vector{uT}
    Dmin::Array{tT,1}
    maxiters::Int64
    step_number::Array{Int64,0}
    initial_extrap::Bool
    length_u::Int64
    tf::tT
    Dtau::tT
end



abstract type IRKAlgorithm{s, initial_extrapolation, mstep} <: OrdinaryDiffEqAlgorithm end
struct IRKGL_seq{s, initial_extrapolation, mstep} <: IRKAlgorithm{s, initial_extrapolation, mstep} end
IRKGL_seq(;s=8, initial_extrapolation=true, mstep=1)=IRKGL_seq{s, initial_extrapolation, mstep}()

function DiffEqBase.__solve(prob::DiffEqBase.AbstractODEProblem{uType,tspanType,isinplace},
    alg::IRKGL_seq{s,initial_extrapolation, mstep}, args...;
    dt=zero(eltype(tspanType)),
    save_on=true,
    adaptive=true,
    maxiters=100,
    kwargs...) where {uType,tspanType,isinplace,s,initial_extrapolation,mstep}


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
    Dtau=convert(tType,dt) # Only used if adaptive=true

#   Memory preallocation (IRKL_Cache)   

    step_number = Array{Int64,0}(undef)
    step_number[] = 0
    init_extrap = initial_extrapolation

    (b, c, a, mu, nu, theta, omega, d) = IRKGLCoefficients(tType,s)
    K = similar(b)
    logK = similar(b)
    Kinv = similar(b)
    Tau = similar(b)
    Tau_ = (1 .+ c)*Dtau
    alpha=similar(theta)
    length_u = length(u0)

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
           
    t0=prob.tspan[1]
    tf=prob.tspan[2] 

    irkgl_cache = IRKGL_Cache(f,p,b,c,a,mu,nu,theta,omega,alpha,d,
                              K,logK,Kinv,Tau,Tau_,
                              U, U_, L, 
                              F,Dmin,maxiters,step_number,
                              init_extrap,length_u,tf,Dtau)
#

    uu = uType[]
    tt = tType[]
                                
    push!(uu,copy(u0))
    push!(tt,t0)
                                                           
    tj = [t0, zero(t0)]
    uj = copy(u0)
    ej= zero(u0)

    if dt==0.
        @warn("Requires a choice of dt>0")
        checks=false
    end

     
    if checks

        cont=true
        error_warn=0
        if !adaptive dt=min(dt,abs(tf-t0)) end
        dts=[dt,dtprev,signdt] 

        if adaptive

            dts[1] = zero(tType)

        
            while cont

                for i in 1:mstep 

                    step_number[]+= 1
                    step_retcode = Main.IRKGLstep_adap!(tj,uj,ej,dts,stats,irkgl_cache)

                    if !step_retcode
                        error_warn=1
                        cont = false
                    end
                    
                    if (tj[1]==tf)  
                        cont=false  
                        break
                    end

                end

                if save_on
                    push!(uu,copy(uj))
                    push!(tt,tj[1])
                end

            end 
        
        else

            while cont

                for i in 1:mstep 

                    step_number[]+= 1
                    step_retcode=Main.IRKGLstep_fixed!(tj,uj,ej,dts,stats,irkgl_cache)

                    if !step_retcode
                        error_warn=1
                        cont = false
                    end
                    
                    if (tj[1]==tf)  
                        cont=false  
                        break
                    end

                end

                if save_on
                    push!(uu,copy(uj))
                    push!(tt,tj[1])
                end



            end 

        end

        stats.naccept=step_number[]
        if error_warn!=0
            @warn("Error during the integration warn=$error_warn")
            sol=DiffEqBase.build_solution(prob,alg,tt,uu,stats=stats,retcode= ReturnCode.Failure)
        
        else

            if !save_on
                push!(uu,copy(uj))
                push!(tt,tj[1])
            end

            sol=DiffEqBase.build_solution(prob,alg,tt,uu,stats=stats,retcode= ReturnCode.Success)
        end

    else

        sol=DiffEqBase.build_solution(prob,alg,tt,uu,stats=stats,retcode= ReturnCode.Failure)

    end

    return(sol)
  
end


