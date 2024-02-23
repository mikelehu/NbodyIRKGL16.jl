@enum retcode begin
    Success
    Failure
end


struct IRKGL_SIMD_Cache{floatT,fType,pType,s_,dim_}
    odef::fType # function defining the ODE system
    p::pType # parameters and so
    b::Vec{s_,floatT}
    c::Vec{s_,floatT}
    mu::VecArray{s_,floatT,2}
    nu::VecArray{s_,floatT,2}
    U::VecArray{s_,floatT,dim_}
    U_::VecArray{s_,floatT,dim_}
    L::VecArray{s_,floatT,dim_}
    F::VecArray{s_,floatT,dim_}
    Dmin::Array{floatT,1}
    maxiters::Int64
    step_number::Array{Int64,0}
    initial_extrap::Array{Int64,0}
    length_u::Int64
    tf::floatT
    step_retcode::retcode
end


abstract type IRKAlgorithm{s, initial_extrap,m, floatType} <: OrdinaryDiffEqAlgorithm end
struct IRKGL_simd{s, initial_extrap,m, floatType} <: IRKAlgorithm{s, initial_extrap,m, floatType} end
IRKGL_simd(;s=8, initial_extrap=1,m=1, floatType=Float64)=IRKGL_simd{s, initial_extrap,m, floatType}()

function DiffEqBase.__solve(prob::DiffEqBase.AbstractODEProblem{uType,tspanType,isinplace},
    alg::IRKGL_simd{s,initial_extrap, m, floatType}, args...;
    dt=0.,
    save_everystep=true,
    adaptive=false,
    maxiters=100,
    kwargs...) where {floatType<: Union{Float32,Float64},uType,tspanType,isinplace,s,initial_extrap, m}

    checks=true

    stats = DiffEqBase.Stats(0)
    stats.nf=0
    stats.nnonliniter=0
    stats.naccept=0

    @unpack f,u0,tspan,p,kwargs=prob
    f= SciMLBase.unwrapped_f(prob.f) 

    tType=eltype(tspanType)
    uiType=eltype(uType)
 
#   Memory preallocation (IRKL_Cache)   

    step_number = Array{Int64,0}(undef)
    step_number[] = 0
    init_interp =  Array{Int64,0}(undef)
    init_interp[] = initial_extrap

    (b_, c_, mu_, nu_) = IRKGLCoefficients(s,dt)
    length_u = length(u0)
    dims = size(u0)
    indices=1:length_u
    step_retcode::retcode=Success

    c = vload(Vec{s,floatType}, c_, 1)
    b = vload(Vec{s,floatType}, b_, 1)
    nu=VecArray{s,floatType,2}(nu_)
    mu=VecArray{s,floatType,2}(mu_)

    zz=zeros(Float64, s, dims...)
    U=VecArray{s,Float64,length(dims)+1}(zz)
    U_=deepcopy(U)
    L=deepcopy(U)
    F=deepcopy(U)

    Dmin=Array{uiType}(undef,length_u)
    Dmin.=zero(uiType)


    dtprev=zero(tType)
    signdt=sign(tspan[2]-tspan[1])    
           
    if signdt==1           # forward integration
        t0=prob.tspan[1]
        tf=prob.tspan[2]   
    else                   # backward integration
        t0=prob.tspan[2]
        tf=prob.tspan[1]   
    end

    dts=[dt,dtprev,signdt]   


    irkgl_cache = IRKGL_SIMD_Cache(f,p,b,c,mu,nu,
                                   U, U_, L, 
                                   F,Dmin,maxiters,step_number,
                                   init_interp,length_u,tf,
                                   step_retcode)

#

    uu = uType[]
    tt = tType[]

    push!(uu,copy(u0))
    push!(tt,t0)
                                                                               
    tj = [t0, zero(t0)]
    uj = copy(u0)
    ej=zero(u0)    
    
                            
    if dt==0.
        println("Error: dt required for fixed timestep.")
        checks=false
    end

    if adaptive==true
        println("Error: only allowed adaptive =false")
        checks=false
    end

    if checks==true

       cont=true
       error_warn=0

       while cont

            for i in 1:m 

                step_number[] += 1
                IRKGLstep_SIMD_fixed!(tj,uj,ej,dts,stats,irkgl_cache)
  
                if (step_retcode==Failure)
                    error_warn=1
                    break
                end
 
                if dts[1]==0 cont=false end
        
            end

            if save_everystep !=false
                push!(uu,uj+ej)
                push!(tt,tj[1]+tj[2])
            end

        end # while

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

