
struct IRKGL_SIMD_Cache{floatT,fType,pType,s_,dim_}
    odef::fType # function defining the ODE system
    p::pType # parameters and so
    b::Vec{s_,floatT}
    c::Vec{s_,floatT}
    a::VecArray{s_,floatT,2}
    mu::VecArray{s_,floatT,2}
    nu::VecArray{s_,floatT,2}
    theta::VecArray{s_,floatT,2}
    omega::VecArray{s_,floatT,2}
    alpha::VecArray{s_,floatT,2}
    d::Vec{s_,floatT}
    s_beta::Vec{s_,floatT}
    K::Vec{s_,floatT}    
    logK::Vec{s_,floatT} 
    Kinv::Vec{s_,floatT}
    Tau::Vec{s_,floatT} 
    Tau_::Vec{s_,floatT}  
    U::VecArray{s_,floatT,dim_}
    U_::VecArray{s_,floatT,dim_}
    L::VecArray{s_,floatT,dim_}
    F::VecArray{s_,floatT,dim_}
    Dmin::Array{floatT,1}
    maxiters::Int64
    step_number::Array{Int64,0}
    initial_extrap::Bool
    length_u::Int64
    tf::floatT
    Dtau::floatT
end


abstract type IRKAlgorithm{s, initial_extrapolation, mstep, floatType} <: OrdinaryDiffEqAlgorithm end
struct IRKGL_simd{s, initial_extrapolation, mstep, floatType} <: IRKAlgorithm{s, initial_extrapolation, mstep, floatType} end
IRKGL_simd(;s=8, initial_extrapolation=true, mstep=1, floatType=Float64)=IRKGL_simd{s, initial_extrapolation, mstep, floatType}()

function DiffEqBase.__solve(prob::DiffEqBase.AbstractODEProblem{uType,tspanType,isinplace},
    alg::IRKGL_simd{s, initial_extrapolation, mstep, floatType}, args...;
    dt=zero(eltype(tspanType)),
    save_on=true,
    adaptive=true,
    maxiters=100,
    kwargs...) where {floatType<: Union{Float32,Float64},uType,tspanType,isinplace,s,initial_extrapolation, mstep}

    checks=true

    stats = DiffEqBase.Stats(0)
    stats.nf=0
    stats.nnonliniter=0
    stats.naccept=0

    @unpack f,u0,tspan,p,kwargs=prob
    f= SciMLBase.unwrapped_f(prob.f) 

    tType=eltype(tspanType)
    uiType=eltype(uType)
    Dtau=convert(tType,dt) # Only used if adaptive=true
 
#   Memory preallocation (IRKL_Cache)   

    step_number = Array{Int64,0}(undef)
    step_number[] = 0
    init_extrap = initial_extrapolation

    (b_, c_, a_, mu_, nu_, theta_, omega_, d_) = IRKGLCoefficients(tType,s)

    length_u = length(u0)
    dims = size(u0)

    b = vload(Vec{s,floatType}, b_, 1)
    c = vload(Vec{s,floatType}, c_, 1)
    a=VecArray{s,floatType,2}(a_)
    nu=VecArray{s,floatType,2}(nu_)
    mu=VecArray{s,floatType,2}(mu_)
    theta=VecArray{s,floatType,2}(theta_)
    omega=VecArray{s,floatType,2}(omega_)
    alpha=deepcopy(theta)
    d = vload(Vec{s,floatType}, d_, 1)
    s_beta= copy(b)

    K = copy(b)
    logK = copy(b)
    Kinv = copy(b)
    Tau = copy(b)
    Tau_ = (1 + c)*Dtau
    zz=zeros(floatType, s, dims...)
    U=VecArray{s,floatType,length(dims)+1}(zz)
    U_=deepcopy(U)
    L=deepcopy(U)
    F=deepcopy(U)

    Dmin=Array{uiType}(undef,length_u)
    Dmin.=zero(uiType)

    dtprev=zero(tType)
    signdt=sign(tspan[2]-tspan[1])    
           
    t0=prob.tspan[1]
    tf=prob.tspan[2] 
 
    irkgl_cache = IRKGL_SIMD_Cache(f,p,b,c,a, mu,nu, theta,omega,alpha,d,s_beta,
                                   K,logK, Kinv, Tau, Tau_,
                                   U, U_, L, 
                                   F,Dmin,maxiters,step_number,
                                   init_extrap,length_u,tf,Dtau)


    uu = uType[]
    tt = tType[]

    push!(uu,copy(u0))
    push!(tt,t0)
                                                                               
    tj = [t0, zero(t0)]
    uj = copy(u0)
    ej=zero(u0)    
                               
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

                    step_number[] += 1
                    step_retcode= Main.IRKGLstep_SIMD_adap!(tj,uj,ej,dts,stats,irkgl_cache)

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

                    step_number[] += 1
                    step_retcode= Main.IRKGLstep_SIMD_fixed!(tj,uj,ej,dts,stats,irkgl_cache)

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

