
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
    length_q::Int64
    tf::tT
    Dtau::tT
end


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
    length_q::Int64
    tf::floatT
    Dtau::floatT
end


abstract type IRKAlgorithm{s, initial_extrapolation, second_order_ode, mstep, simd, floatType} <: OrdinaryDiffEq.OrdinaryDiffEqAlgorithm end
struct nbirkgl16{s, initial_extrapolation, second_order_ode, mstep, simd, floatType} <: IRKAlgorithm{s, initial_extrapolation, second_order_ode, mstep, simd, floatType} end
nbirkgl16(;s=8, initial_extrapolation=true, second_order_ode=true, mstep=1, simd=true, floatType=Float64)=nbirkgl16{s, initial_extrapolation, second_order_ode, mstep, simd,floatType}()

function DiffEqBase.__solve(prob::DiffEqBase.AbstractODEProblem{uType,tspanType,isinplace},
    alg::nbirkgl16{s,initial_extrapolation, second_order_ode, mstep, simd, floatType}, args...;
    dt=zero(eltype(tspanType)),
    save_everystep=true,
    adaptive=true,
    maxiters=100,
    kwargs...) where {floatType<:Number, uType,tspanType,isinplace,s,initial_extrapolation,second_order_ode,mstep, simd}


    checks=true

    stats = DiffEqBase.Stats(0)
    stats.nf=0
    stats.nfpiter=0
    stats.naccept=0
    
    @unpack f,u0,tspan,p,kwargs=prob
    f= SciMLBase.unwrapped_f(prob.f) 

    t0=prob.tspan[1]
    tf=prob.tspan[2] 

    length_u = length(u0)
    length_q=div(length_u,2)

    step_fun::Function=empty

    if simd

        if second_order_ode
            if adaptive
                step_fun=Main.Main.IRKNGLstep_SIMD_adap_simpl!
            else
                step_fun=Main.Main.IRKNGLstep_SIMD_fixed_simpl!
            end
    
        else
            if adaptive
                step_fun=Main.IRKGLstep_SIMD_adap!
            else
                step_fun=Main.IRKGLstep_SIMD_fixed!
            end
    
        end

    else

        if second_order_ode
            if adaptive
                step_fun=Main.IRKNGLstep_adap_simpl!
            else
                step_fun=Main.IRKNGLstep_fixed_simpl!
            end

        else
            if adaptive
                step_fun=Main.IRKGLstep_adap!
            else
                step_fun=Main.IRKGLstep_fixed!
            end

        end

    end

    
    uiType=eltype(uType)
    if uiType <:Complex 
        uiType=real(uiType)
    end

    tType=eltype(tspanType)
    Dtau=convert(tType,dt) # Only used if adaptive=true

    if dt==0.
        @warn("Requires a choice of dt>0")
        checks=false
    end

    if simd
        if !(uiType <: Union{Float32,Float64}) || !(tType <:Union{Float32,Float64}) || !(floatType <:Union{Float32,Float64})
            @warn("SIMD vectorization can be used only with Float64 or Float32. \n Please try again with simd=false option")
            checks=false
        end
    end
     
    if checks

        #   Memory preallocation (IRKL_Cache)   

        Dmin=Array{uiType}(undef,length_u)
        Dmin.=zero(uiType)

        step_number = Array{Int64,0}(undef)
        step_number[] = 0
        init_extrap = initial_extrapolation

        if simd

            (b_, c_, a_, mu_, nu_, theta_, omega_, d_) = IRKGLCoefficients(tType,s)
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
           
            irkgl_cache = IRKGL_SIMD_Cache(f,p,b,c,a, mu,nu, theta,omega,alpha,d,s_beta,
                                        K,logK, Kinv, Tau, Tau_,
                                        U, U_, L, 
                                        F,Dmin,maxiters,step_number,
                                        init_extrap,length_u, length_q, tf,Dtau)
    

        else

            (b, c, a, mu, nu, theta, omega, d) = IRKGLCoefficients(tType,s)
            K = similar(b)
            logK = similar(b)
            Kinv = similar(b)
            Tau = similar(b)
            Tau_ = (1 .+ c)*Dtau
            alpha=similar(theta)

            U=Array{uType}(undef, s)
            for i in 1:s
                U[i]=zero(u0)
            end

            U_=deepcopy(U)
            L=deepcopy(U)
            F=deepcopy(U)

    
            irkgl_cache = IRKGL_Cache(f,p,b,c,a,mu,nu,theta,omega,alpha,d,
                                    K,logK,Kinv,Tau,Tau_,
                                    U, U_, L, 
                                    F,Dmin,maxiters,step_number,
                                    init_extrap,length_u, length_q, tf,Dtau)
        
        end

        uu = uType[]
        tt = tType[]
                                    
        push!(uu,copy(u0))
        push!(tt,t0)
                                                            
        tj = [t0, zero(t0)]
        uj = copy(u0)
        ej= zero(u0)

        dtprev=zero(tType)
        signdt=sign(tspan[2]-tspan[1]) 

        if !adaptive dt=min(abs(dt),abs(tf-t0)) end
        dts=[dt,dtprev,signdt] 

        if adaptive dts[1] = zero(tType) end

        cont=true
        error_warn=0


        while cont

                for i in 1:mstep 

                    step_number[]+= 1
                    step_retcode = step_fun(tj,uj,ej,dts,stats,irkgl_cache)

                    if !step_retcode
                        error_warn=1
                        cont = false
                    end
                    
                    if (tj[1]==tf)  
                        cont=false  
                        break
                    end

                end

                if save_everystep
                    push!(uu,copy(uj))
                    push!(tt,tj[1])
                end

        end 
        
        stats.naccept=step_number[]

        if error_warn!=0
            @warn("Error during the integration warn=$error_warn")
            sol=DiffEqBase.build_solution(prob,alg,tt,uu,stats=stats,retcode= ReturnCode.Failure)
        
        else

            if !save_everystep
                push!(uu,copy(uj))
                push!(tt,tj[1])
            end

            sol=DiffEqBase.build_solution(prob,alg,tt,uu,stats=stats,retcode= ReturnCode.Success)
        end

    else

        sol=DiffEqBase.build_solution(prob,alg,[],[],stats=stats,retcode= ReturnCode.Failure)

    end

    return(sol)
  
end


