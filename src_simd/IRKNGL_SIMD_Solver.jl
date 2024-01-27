

include("VecArray_def.jl")


struct IRKNGL_SIMD_Cache{floatType,fType,pType,s_,dim,dim_}
    odef::fType # function defining the ODE system
    p::pType # parameters and so
    b::Vec{s_,floatType}
    c::Vec{s_,floatType}
    mu::VecArray{s_,floatType,2}
    nu::VecArray{s_,floatType,2}
    theta::VecArray{s_,floatType,2}
    omega::VecArray{s_,floatType,2}
    alpha::VecArray{s_,floatType,2}
    g::Vec{s_,floatType}
    d::Vec{s_,floatType}
    s_beta::Vec{s_,floatType}
    U::VecArray{s_,floatType,dim_}
    U_::VecArray{s_,floatType,dim_}
    L::VecArray{s_,floatType,dim_}
    dU::VecArray{s_,floatType,dim_}
    F::VecArray{s_,floatType,dim_}
    Dmin::Array{floatType,dim}
    maxiters::Int64
    step_number::Array{Int64,0}
    tau::Array{floatType,0}
    initial_interp::Array{Int64,0}
    length_u::Int64
    length_q::Int64
    nrmdigits::Array{floatType, 0}
    E_weights::Array{floatType,dim}
end


function Rdigits(x::Real,r::Real)
    mx=r*x
    mxx=mx+x
    return mxx-mx
end

function Rdigits(x::Vec,r::Real)
    mx=r*x
    mxx=mx+x
    return mxx-mx
end


abstract type IRKAlgorithm{s,initial_interp, dim,floatType,m,myoutputs,nrmbits} <: OrdinaryDiffEqAlgorithm end
struct IRKNGL_simd{s, initial_interp, dim,floatType,m,myoutputs,nrmbits} <: IRKAlgorithm{s, initial_interp, dim,floatType,m,myoutputs,nrmbits} end
IRKNGL_simd(;s=8, initial_interp=1, dim=1,floatType=Float64,m=1,myoutputs=false,nrmbits=0)=IRKNGL_simd{s, initial_interp, dim,floatType,m,myoutputs,nrmbits}()

function DiffEqBase.__solve(prob::DiffEqBase.AbstractODEProblem{uType,tspanType,isinplace},
    alg::IRKNGL_simd{s,initial_interp, dim,floatType, m,myoutputs,nrmbits}, args...;
    dt,
    saveat=eltype(tspanType)[],
    gamma=eltype(uType)[],
    save_everystep=true,
    adaptive=false,
    maxiters=100,
    kwargs...) where {floatType<: Union{Float32,Float64},uType,tspanType,isinplace,dim,s,m,initial_interp,myoutputs,nrmbits}

    trace=false
   
    #destats = DiffEqBase.DEStats(0)
    stats = DiffEqBase.Stats(0)
    
    @unpack f,u0,tspan,p,kwargs=prob
    f= SciMLBase.unwrapped_f(prob.f) 

    tType=eltype(tspanType)

    utype = Vector{floatType}
#    ttype = floatType

    step_number = Array{Int64,0}(undef)
    step_number[] = 0
    init_interp =  Array{Int64,0}(undef)
    init_interp[] = initial_interp

#    dts = Array{tType}(undef, 1)
    uiType=eltype(u0)

    dtprev=zero(tType)
    signdt=sign(tspan[2]-tspan[1]) 
               
    if signdt==1 
        t0=prob.tspan[1]
        tf=prob.tspan[2]   # forward
     else
        t0=prob.tspan[2]
        tf=prob.tspan[1]   # backward
     end

    dts=[dt,dtprev,signdt]

    (b_, c_, mu_, nu_, theta_, omega_, g_, d_) = IRKGLCoefficients(s,dt)
    length_u = length(u0)
    length_q=div(length_u,2)
    dims = size(u0)
    indices=1:length_u
    indices2 = (length_q+1):length_u

    c = vload(Vec{s,floatType}, c_, 1)
    b = vload(Vec{s,floatType}, b_, 1)
    nu=VecArray{s,floatType,2}(nu_)
    mu=VecArray{s,floatType,2}(mu_)
    theta=VecArray{s,floatType,2}(theta_)
    omega=VecArray{s,floatType,2}(omega_)
    alpha=deepcopy(theta)
    g = vload(Vec{s,floatType}, g_, 1)
    d = vload(Vec{s,floatType}, d_, 1)
    s_beta= deepcopy(g)


    uu = uType[]
    tt = tType[]

    zz=zeros(Float64, s, dims...)
    U=VecArray{s,Float64,length(dims)+1}(zz)
    U_=deepcopy(U)
    L=deepcopy(U)
    dU=deepcopy(U)
    F=deepcopy(U)

    Dmin=Array{uiType}(undef,length_q)
    Dmin.=zero(uiType)

    E_weights=Array{uiType}(undef,length_u)
    if isempty(gamma)
        E_weights.=one(uiType)
     else
        for k in indices
             E_weights[k]=gamma[k]
        end
     end


    tau=Array{uiType,0}(undef)
    tau[]=zero(uiType)


    nrmdig = Array{uiType, 0}(undef)
    if (nrmbits > 0)
        nrmdig[] = uiType(2^nrmbits)
    else
        nrmdig[] = zero(uiType)
    end

    irkngl_cache = IRKNGL_SIMD_Cache(f,p,b,c,mu,nu,theta,omega,alpha,g,d,s_beta,
                                     U, U_, L, dU,
                                     F,Dmin,maxiters,step_number,tau, 
                                     init_interp,length_u,length_q,nrmdig, E_weights)

    iters = Float64[]
    push!(uu,copy(u0))
    push!(tt,t0)
    push!(iters,0.)
                                                  
    tj = [t0, zero(t0)]
    uj = copy(u0)
    ej=zero(u0)
    
    tj_=similar(tj)
    uj_=similar(uj)
    ej_=similar(ej)
    L_=deepcopy(L)
    dts_=similar(dts)

    tstops=tType[]

    if isempty(saveat)     

        tstops=tType[Inf]  
#        save_everystep=true

    else

        if saveat isa Number 

            if t0<tf
                tstops=vcat(Inf, abs.(reverse(t0:saveat:tf)))
            else
                tstops=vcat(Inf, abs.(tf:-saveat:t0))
            end
            
        else
            tstops=vcat(Inf,abs.(reverse(saveat)))
        end

#       m=1
        save_everystep=false

        if tstops[end]==t0 pop!(tstops) end 

    end

    tout= Array{uiType,0}(undef)
    tout[]=pop!(tstops)
    save_step=false  

    if adaptive

        cont=true

        while cont

            tit=0
            j=0

            trace ? println("Urratsa j=",step_number[],",tj=", tj[1]) : nothing

            for i in 1:m
            
                j+=1
                step_number[] += 1

                if step_number[]==1

                    (status,j_iter) = IRKNGLstep_SIMD_fixed!(tj,tf, uj,ej,prob,dts,irkngl_cache)

                    E2=zero(uiType)
                    for k in indices2
                        Fk=getindex_(F,k) 
                        Ek=dt*sum(g*Fk)      # Ek = dt*dot(g, Fk)
                        wEk=E_weights[k]*Ek
                        E2=muladd(Ek,wEk,E2)
                    end
                    tau[]=E2^(1/(2s-2)) 

                    trace ? println("tau=", tau[], ",E2=", E2) : nothing

                 else
                    (status,j_iter) = IRKNGLstep_SIMD_adap!(tj,tf,uj,ej,prob,dts,irkngl_cache)
                 end

                 if (status=="Failure")
                    println("Fail")
                    sol=DiffEqBase.build_solution(prob,alg,tt,uu, retcode= ReturnCode.Failure)
                    if (myoutputs==true)
                        return(sol,iters, step_number[])
                    else
                        return(sol)
                    end
                end

                tit+=j_iter

                if tout[]<tj[1]+tj[2]

                    tj_.=tj
                    uj_.=uj
                    ej_.=ej
                    L_.data.=L.data
                    init_interp_=init_interp

                    dts_[1]=abs(tj[1]+tj[2]-tout[])
                    dts_[2]=zero(dts[2])
                    dts_[3]=-dts[3]

                    init_interp=0
                    (status_,j_iter_) = IRKNGLstep_SIMD_fixed!(tj_,tf, uj_,ej_,prob,dts_,irkngl_cache)
                    
                    if (status=="Failure")
                        println("Fail- Computing tout")
                        sol=DiffEqBase.build_solution(prob,alg,tt,uu, retcode= ReturnCode.Failure)
                        if (myoutputs==true)
                          return(sol,iters, step_number[])
                        else
                          return(sol)
                        end
                    end

                    L.data.=L_.data
                    init_interp=init_interp_

                    push!(uu,uj_+ej_)
                    push!(tt,tj_[1]+tj_[2])

                    tout[]=pop!(tstops)

                elseif tout[]==tj[1]+tj[2]
                    save_step=true
                    tout[]=pop!(tstops)

                end

                if (dts[1]==0)
                    cont=false
                    break
                end

            end  # for 

            if save_step==true || save_everystep !=false || (cont==false)
                push!(uu,uj+ej)
                push!(tt,tj[1]+tj[2])
                push!(iters, tit/j)
                save_step=false
            end

        end #end while

#        sol=DiffEqBase.build_solution(prob,alg,tt,uu,destats=destats,retcode= ReturnCode.Success)
        sol=DiffEqBase.build_solution(prob,alg,tt,uu,stats=stats,retcode= ReturnCode.Success)

        if (myoutputs==true)
            return(sol,iters, step_number[])
        else
            return(sol)
        end


    else # adaptive=false
                                    

        cont=true
                               
        while cont

            tit=0
            j=0

            for i in 1:m
            
                 j+=1

                  step_number[] += 1
                 (status,j_iter) = IRKNGLstep_SIMD_fixed!(tj,tf,uj,ej,prob,dts,irkngl_cache)
  
                 if (status=="Failure")
                    println("Fail")
                    sol=DiffEqBase.build_solution(prob,alg,tt,uu,retcode= ReturnCode.Failure) 
                    if (myoutputs==true)
                        return(sol,iters, step_number[])
                    else
                        return(sol)
                    end
                 end
  
                 tit+= j_iter

                 if tout[]<tj[1]+tj[2]

                    tj_.=tj
                    uj_.=uj
                    ej_.=ej
                    L_.data.=L.data
                    init_interp_=init_interp

                    dts_[1]=abs(tj[1]+tj[2]-tout[])
                    dts_[2]=zero(dts[2])
                    dts_[3]=-dts[3]

                    init_interp=0
                    (status_,j_iter_) = IRKNGLstep_SIMD_fixed!(tj_,tf, uj_,ej_,prob,dts_,irkngl_cache)
                    
                    if (status=="Failure")
                        println("Fail- Computing tout")
                        sol=DiffEqBase.build_solution(prob,alg,tt,uu, retcode= ReturnCode.Failure)
                        if (myoutputs==true)
                          return(sol,iters, step_number[])
                        else
                          return(sol)
                        end
                    end

                    L.data.=L_.data
                    init_interp=init_interp_

                    push!(uu,uj_+ej_)
                    push!(tt,tj_[1]+tj_[2])

                    tit+=j_iter_-j_iter
                    push!(iters, tit/j)
                    j=0
                    tit=0

                    tout[]=pop!(tstops)

                elseif tout[]==tj[1]+tj[2]
                    save_step=true
                    tout[]=pop!(tstops)

                end
  
                if (dts[1]==0)
                     cont=false
                     break
                end
            
            end
  
            if save_step==true || save_everystep !=false || (cont==false)
                push!(uu,uj+ej)
                push!(tt,tj[1]+tj[2])
                push!(iters, tit/j)
                save_step=false
            end

        end # while

        dts[1]=dt
        dts[2]=zero(tType)

#        sol=DiffEqBase.build_solution(prob,alg,tt,uu,destats=destats,retcode= ReturnCode.Success)
        sol=DiffEqBase.build_solution(prob,alg,tt,uu,stats=stats,retcode= ReturnCode.Success)
  
        if (myoutputs==true)
            return(sol,iters,step_number[])
        else
            return(sol)
        end
 
    end # adaptive-else

end

