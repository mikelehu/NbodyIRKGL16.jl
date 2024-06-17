# NbodyIRKGL16.jl

## **Few-body integrator with time-reversible adaptivity in Julia**

Time-reversible or symplectic methods are very useful for long-term integration of planetary systems with a constant step size. However, in the case of n-body problems involving close encounters, adaptive time-stepping is required. We present a SIMD-vectorized integrator for few-body problems which incorporates a robust time-reversible adaptivity mechanism that makes it highly performant for the long-term numerical integration of problems with close-encounters.


## Description


Time-reversible or symplectic methods are very useful for long-term integration of planetary systems with a constant step size. The 16th order method implemented in our package SciML/IRKGaussLegendre.jl is both time-reversible and symplectic, and includes a SIMD-vectorized version IRKGL16-SIMD (presented in Juliacon2020) that is particularly efficient for problems without close encounters. However, numerical precision degrades greatly during close encounters if constant step size is employed. Close encounters are common, for instance, in simulations of the orbital evolution of the Solar System when some asteroids or comets are included in the model. Hence, a robust numerical integrator of few-body problems must include some adaptive mechanism to deal with close encounters. Unfortunately, with conventional time-step adaptivity strategies the advantages of symplectic methods is lost. 

We present an integrator for few-body problems, which based on IRKGL16-SIMD, that incorporates a robust time-reversible adaptivity mechanism that makes it highly performant for the long-term numerical integration of problems with close-encounters. Fortunately, time-stepping strategies that preserve the reversible structure of the n-body problem are possible, for instance, by applying reversible stepsize strategies introduced in [2], which require an explicitly given expression of the stepsize as a function of the positions and velocities. In our code, we make use of an appropriate step-size function for the n-problem proposed in our previous work [1].

We show some numerical experiments with few-body problems that demonstrates the efficiency and robustness of our code.


## References


- [1] Global Time-Renormalization of the Gravitational N-body Problem, M. Antoñana, P. Chartier, J. Makazaga and A. Murua. SIAM Journal on Applied Dynamical System (2020). https://doi.org/10.1137/20M1314719.

- [2] Reversible Long-Term Integration with Variable Stepsizes,  E.Hairer and  D. Stoffer.
SIAM Journal on Scientific Computing (1997).
https://doi.org/10.1137/S1064827595285494

## Repository

-https://github.com/mikelehu/NbodyIRKGL16.jl



