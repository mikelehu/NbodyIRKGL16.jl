# NbodyIRKGL16.jl :  **few-body integrator with time-reversible adaptivity in Julia**



## Description


We present an integrator for few-body problems, which based on IRKGL16-SIMD, that incorporates a robust time-reversible adaptivity mechanism that makes it highly performant for the long-term numerical integration of problems with close-encounters. 



## Installation

```julia
PATH_SRC_SIMD="../../src/simd/"
include(string(PATH_SRC_SIMD,"IRKGL_SIMD.jl"))
using .IRKGL_SIMD   
```

## Solver options


### Available common arguments


- dt:
    - if adaptive=false, dt is the stepsize for the integration
    - if adaptive=true, dt is a constant to specify the tolerance (Dtau=dt)

- save_on: denotes whether intermediate solutions are saved (default is true)
- adaptive =true (adaptive timestepping); =false (fixed timestepping)
- maxiters: maximum number of iterations before stopping


### No-common arguments


- initial_extrapolation: initialization method for stages:

        - =false  simplest initialization
        - =true (default) interpolating from the stage values of previous step

- mstep: output saved at every 'mstep' steps (default mstep=1)



## Return Codes

The solution types have a retcode field which returns a symbol signifying the error state of the solution. The retcodes are as follows:

- ReturnCode.Success: The integration completed without errors.
- ReturnCode.Failure: General uncategorized failures or errors.


## Example: 15-body model of the Solar System


We consider the 15-body model of the Solar Systen that includes: the Sun, all eight planets of the Solar System, Pluto and the five main bodies of the asteroid belt.

We consider the initial values at Julian day (TDB) 2440400.5 (the 28th of June of 1969), obtained form the DE430 ephemerides. 

We run a numerical integration back in time for $t=2e7$ days where a close approach between Ceres and Bamberga asteroids occurs.

### Step 1: Defining  the problem

To solve this numerically, we define a problem type by giving it the equation, the initial
condition, and the timespan to solve over:

```julia
using NbodyIRKGL16
using Plots, LinearAlgebra

PATH_ODES="../../ODEProblems/"
include(string(PATH_ODES,"Initial15Body.jl"))
include(string(PATH_ODES,"Nbody.jl"));
```

```julia
function NbodyODE_fstep!(F,u,Gm,t,part)

   N = length(Gm)

   if part==1  # Evaluate dq/dt

      for i in 1:3, j in 1:N
         F[i,j,1] = u[i,j,2]
      end

      sinv=zero(eltype(typeof(u)))
   
   else        # Evaluate dv/dt

      kappa=3^2   

      A = zero(eltype(u))
      B = zero(eltype(u))
      C = zero(eltype(u))

      for i in 1:N
         for k in 1:3
            F[k, i, 2] = 0
         end
      end

      for i in 1:N
         xi = u[1,i,1]
         yi = u[2,i,1]
         zi = u[3,i,1]
         vxi = u[1,i,2]
         vyi = u[2,i,2]
         vzi = u[3,i,2]
         Gmi = Gm[i]
         for j in i+1:N
            xij = xi - u[1,j,1]
            yij = yi - u[2,j,1]
            zij = zi - u[3,j,1]
            vxij = vxi - u[1,j,2]
            vyij = vyi - u[2,j,2]
            vzij = vzi - u[3,j,2]
            Gmj = Gm[j]
            invnorm2qij =1/(xij*xij+yij*yij+zij*zij)
            invnormqij = sqrt(invnorm2qij)
            auxij = invnorm2qij * invnormqij
            Gmjauxij = Gmj*auxij
            F[1,i,2] -= Gmjauxij*xij
            F[2,i,2] -= Gmjauxij*yij
            F[3,i,2] -= Gmjauxij*zij
            Gmiauxij = Gmi*auxij
            F[1,j,2] += Gmiauxij*xij
            F[2,j,2] += Gmiauxij*yij
            F[3,j,2] += Gmiauxij*zij
            norm2vij = vxij*vxij+vyij*vyij+vzij*vzij
            A += (norm2vij*invnorm2qij)^2
            B += (Gmi+Gmj)*invnorm2qij
            C += invnorm2qij 
         end
      end

      sinv = (kappa*A  + B^2 * C)^(1/4)    # 0.25
   
   end

   return sinv

end
```

```julia

fltype=Float64
u0, Gm, bodylist = Initial15Body(fltype)  # defined in "Initial15Body.jl" file.
N = length(Gm)
show(bodylist)

t0=fltype(0.)
tF=fltype(-2e4)  
tspan= (t0,tF)
prob = ODEProblem(NbodyODE_fstep!, u0,tspan , Gm);

```

### Step 2: Solving the problem


After defining a problem, you solve it using solve

```julia
Dtau=fltype(1.8)
alg= IRKNGL_simd(initial_extrapolation=true)
sol=solve(prob, alg, adaptive=true, dt = Dtau)
sol.retcode
```

### Step 3: Analyzing the solution


#### Orbits of planets and asteroids


```julia
pl1 = plot(title="Inner Solar System", 
            xlabel="x", ylabel="y",  aspect_ratio=1)

for j = 2:5
    x  = [u[1,j,1] for u in sol.u]
    y  = [u[2,j,1] for u in sol.u] 
    pl1 = plot!(x,y, label=bodylist[j]) 
end 


pl2 = plot(title="Outer Solar System", 
            xlabel="x", ylabel="y",  aspect_ratio=1)

for j = 6:10
    x  = [u[1,j,1] for u in sol.u]
    y  = [u[2,j,1] for u in sol.u] 
    pl2 = plot!(x,y, label=bodylist[j]) 
end 


pl3 = plot(title="Big asteroids", 
            xlabel="x", ylabel="y",  aspect_ratio=1)

for j = 11:15
    x  = [u[1,j,1] for u in sol.u]
    y  = [u[2,j,1] for u in sol.u] 
    pl3 = plot!(x,y, 
        label="")
end 


plot(pl1,pl2,pl3, layout=(1,3), size=(1200,300))
```
![15-body Solar System](/Examples/BodyOrbits.png)



#### Error in Energy

```julia
function NbodyEnergy(u,Gm)

     N = length(Gm)
     zerouel = zero(eltype(u))
     T = zerouel
     U = zerouel
     for i in 1:N
        qi = u[2,:,i]
        vi = u[1,:,i]
        Gmi = Gm[i]
        T += Gmi*(vi[1]*vi[1]+vi[2]*vi[2]+vi[3]*vi[3])
        for j in (i+1):N
           qj = u[2,:,j]  
           Gmj = Gm[j]
           qij = qi - qj
           U -= Gmi*Gmj/norm(qij)
        end
     end

    1/2*T + U

end
```


```julia
yrange=(1e-18,1e-14)

E0=NbodyEnergy(BigFloat.(u0), BigFloat.(Gm))
ΔE = map(x->NbodyEnergy(BigFloat.(x),BigFloat.(Gm)), sol.u)./E0.-1;


pl1=plot(title="Error in Energy", xlabel="t", ylabel="log10(ΔE/E0)",
        titlefontsize=18,
        xtickfont = font(12),
        ytickfont = font(12),
        xticks=([-2e4,-1.5e4,-1e4,-5e3, 0],["-2e4","-1.5e4","-1e4","-5e3", "0"]),
        yscale=:log10,
        legend=:topright, ylims=yrange)

plot!(pl1,sol.t[2:end-1],abs.(ΔE[2:end-1]), color=:red, label="")
```

![Error in energy](/Examples/EnergyError.png)


#### Close encounters between asteroids

```julia
n=length(sol.t)
dist=Array{Vector{Float64}}(undef,5,5)
for i in 1:5
    for j in 1:5
        dist[i,j]=zeros(n)
    end
end

for i in 1:5
    ix=10+i        # indices for asteroids: 11:15
    A=[u[:,ix,1]  for u in sol.u]
    for j in i+1:5
        jx=10+j   # indices for asteroids: 11:15
        B=[u[:,jx,1]  for u in sol.u]
        dist[i,j]=map(x->norm(x), A-B)
    end
end
```

```julia
pl1 = plot(title="Pairwise distances between asteroids", 
           titlefontsize=18,
           xtickfont = font(12),
           ytickfont = font(12),
           yscale=:log10,
           legend=:bottomright,
                   xticks=([-2e4,-1.5e4,-1e4,-5e3, 0],["-2e4","-1.5e4","-1e4","-5e3", "0"]),
           xlabel="t", ylabel="log10(dist)")


for i in 1:5
    for j in i+1:5
        if (i==1 && j==5)
             plot!(pl1, sol.t,dist[i,j], label="Ceres-Bamberga")
        else
            plot!(pl1, sol.t,dist[i,j], label="")
        end
    end
end



pl2= scatter(title="Step size",
             titlefontsize=18,
             xtickfont = font(12),
             ytickfont = font(12),
             xticks=([-2e4,-1.5e4,-1e4,-5e3, 0],["-2e4","-1.5e4","-1e4","-5e3", "0"]),
             xlabel="t", ylabel="h",
             label="",
             sol.t[1:end-2],
             abs.(sol.t[2:end-1].-sol.t[1:end-2]))


plot(pl1,pl2, layout=(1,2), size=(900,300))
```


![Close encounter](/Examples/CloseEncounter.png)


## References


- [1] Global Time-Renormalization of the Gravitational N-body Problem, M. Antoñana, P. Chartier, J. Makazaga and A. Murua. SIAM Journal on Applied Dynamical System (2020). https://doi.org/10.1137/20M1314719.

- [2] Reversible Long-Term Integration with Variable Stepsizes,  E.Hairer and  D. Stoffer.
SIAM Journal on Scientific Computing (1997).
https://doi.org/10.1137/S1064827595285494

## Repository

-https://github.com/mikelehu/NbodyIRKGL16.jl



