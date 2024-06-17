#
# NbodyEnergy(u,Gm)
# NbodyBarycenter(u,Gm)
# Compute the right side of the ODE system of the N-body problem:
#     - NbodyODE!(F,u,Gm,t)
#     - NbodyODE!(F,u,Gm,t,part)
#
# Same as NbodyODE!() but it addtion evaluate and return the step-size expression: 
#     - NbodyODE_fstep!(F,u,Gm,t) 
#     - NbodyODE_fstep!(F,u,Gm,t,part) 
#


function NbodyEnergy(u,Gm)

     N = length(Gm)
     zerouel = zero(eltype(u))
     T = zerouel
     U = zerouel

     for i in 1:N
        qi = u[:,i,1]
        vi = u[:,i,2]
        Gmi = Gm[i]
        T += Gmi*(vi[1]*vi[1]+vi[2]*vi[2]+vi[3]*vi[3])
        for j in (i+1):N
           qj = u[:,j,1]
           Gmj = Gm[j]
           qij = qi - qj
           U -= Gmi*Gmj/norm(qij)
        end
     end

    1/2*T + U

end



function NbodyBarycenter(u,Gm)

     N = length(Gm)
     dim = size(u,1)
     A = zeros(dim)
     B = zeros(dim)

     for i in 1:N
        qi = u[:,i,1]
        vi = u[:,i,2]
        Gmi = Gm[i]
        A += Gmi*qi
        B += Gmi*vi
     end

     return A, B

end



function NbodyODE!(F,u,Gm,t)

     N = length(Gm)

     for i in 1:N
        for k in 1:3
            F[k, i, 2] = 0
        end
     end

     for i in 1:N
        xi = u[1,i,1]
        yi = u[2,i,1]
        zi = u[3,i,1]
        Gmi = Gm[i]
        for j in i+1:N
            xij = xi - u[1,j,1]
            yij = yi - u[2,j,1]
            zij = zi - u[3,j,1]
            Gmj = Gm[j]
            dotij = (xij*xij+yij*yij+zij*zij)
            auxij = 1/(sqrt(dotij)*dotij)
            Gmjauxij = Gmj*auxij
            F[1,i,2] -= Gmjauxij*xij
            F[2,i,2] -= Gmjauxij*yij
            F[3,i,2] -= Gmjauxij*zij
            Gmiauxij = Gmi*auxij
            F[1,j,2] += Gmiauxij*xij
            F[2,j,2] += Gmiauxij*yij
            F[3,j,2] += Gmiauxij*zij
        end
     end

     for i in 1:3, j in 1:N
        F[i,j,1] = u[i,j,2]
     end

    return nothing
end


function NbodyODE!(F,u,Gm,t,part)
    
   
   N = length(Gm)

  if part==1    # Evaluate dq/dt

       for i in 1:3, j in 1:N
          F[i,j,1] = u[i,j,2]
       end
  else         # Evaluate dv/dt
       for i in 1:N
          for k in 1:3
              F[k, i, 2] = 0
          end
       end

       for i in 1:N
          xi = u[1,i,1]
          yi = u[2,i,1]
          zi = u[3,i,1]
          Gmi = Gm[i]
          for j in i+1:N
                 xij = xi - u[1,j,1]
                 yij = yi - u[2,j,1]
                 zij = zi - u[3,j,1]
                 Gmj = Gm[j]
                 dotij = (xij*xij+yij*yij+zij*zij)
                 auxij = 1/(sqrt(dotij)*dotij)
                 Gmjauxij = Gmj*auxij
                 F[1,i,2] -= Gmjauxij*xij
                 F[2,i,2] -= Gmjauxij*yij
                 F[3,i,2] -= Gmjauxij*zij
                 Gmiauxij = Gmi*auxij
                 F[1,j,2] += Gmiauxij*xij
                 F[2,j,2] += Gmiauxij*yij
                 F[3,j,2] += Gmiauxij*zij
          end
       end
  end 

  return nothing

end



function NbodyODE_fstep!(F,u,Gm,t)

   kappa=3^2   # kappa=alpha^p, alpha=3, p=2
   N = length(Gm)
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

   for i in 1:3, j in 1:N
      F[i,j,1] = u[i,j,2]
   end

   sinv = (kappa*A  + B^2 * C)^(1/4)    # 0.25
   return sinv

end


function NbodyODE_fstep!(F,u,Gm,t,part)

   N = length(Gm)

   if part==1  # Evaluate dq/dt

      for i in 1:3, j in 1:N
         F[i,j,1] = u[i,j,2]
      end

      sinv=zero(eltype(typeof(u)))
   
   else        # Evaluate dv/dt

      kappa=3^2   # kappa=alpha^p, alpha=3, p=2

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



