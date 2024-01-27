



function NbodyEnergy3REL(u,Gm)

     N = length(Gm)

     norm12=sqrt(u[1,1,1]^2+u[2,1,1]^2)
     norm23=sqrt(u[1,2,1]^2+u[2,2,1]^2)
     norm31=sqrt(u[1,3,1]^2+u[2,3,1]^2)

     T=Gm[1]*(u[1,1,2]^2+u[2,1,2]^2)+Gm[2]*(u[1,2,2]^2+u[2,2,2]^2)+Gm[3]*(u[1,3,2]^2+u[2,3,2]^2)
     U=Gm[1]*Gm[2]/norm12+Gm[2]*Gm[3]/norm23+Gm[3]*Gm[1]/norm31

     return T/2-U


end


function NbodyODE3REL!(du,u,Gm,t)
     
      N = length(Gm)
      
     
      q12_x = u[1,1,1]
      q12_y = u[2,1,1]
      q23_x = u[1,2,1]
      q23_y = u[2,2,1]
      q31_x = u[1,3,1]
      q31_y = u[2,3,1]

      dot12=(q12_x*q12_x+q12_y*q12_y)^-1.5
      dot23=(q23_x*q23_x+q23_y*q23_y)^-1.5
      dot31=(q31_x*q31_x+q31_y*q31_y)^-1.5
      
      aux12_x = q12_x*dot12
      aux12_y = q12_y*dot12
      aux23_x = q23_x*dot23
      aux23_y = q23_y*dot23
      aux31_x = q31_x*dot31
      aux31_y = q31_y*dot31

      dv1_x = Gm[2]*aux12_x - Gm[3]*aux31_x
      dv1_y = Gm[2]*aux12_y - Gm[3]*aux31_y
      dv2_x = Gm[3]*aux23_x - Gm[1]*aux12_x 
      dv2_y = Gm[3]*aux23_y - Gm[1]*aux12_y 
      dv3_x = Gm[1]*aux31_x - Gm[2]*aux23_x 
      dv3_y = Gm[1]*aux31_y - Gm[2]*aux23_y 

      du[1,1,1]=u[1,2,2]-u[1,1,2]
      du[2,1,1]=u[2,2,2]-u[2,1,2]
      du[1,2,1]=u[1,3,2]-u[1,2,2]
      du[2,2,1]=u[2,3,2]-u[2,2,2]
      du[1,3,1]=u[1,1,2]-u[1,3,2]
      du[2,3,1]=u[2,1,2]-u[2,3,2]

      du[1,1,2] = dv1_x
      du[2,1,2] = dv1_y
      du[1,2,2] = dv2_x
      du[2,2,2] = dv2_y
      du[1,3,2] = dv3_x 
      du[2,3,2] = dv3_y
        

    return nothing

end




function NbodyODE3REL!(du,u,Gm,t,part)
    
   
     N = length(Gm)

    if part==1    # Evaluate dq/dt

      du[1,1,1]=u[1,2,2]-u[1,1,2]
      du[2,1,1]=u[2,2,2]-u[2,1,2]
      du[1,2,1]=u[1,3,2]-u[1,2,2]
      du[2,2,1]=u[2,3,2]-u[2,2,2]
      du[1,3,1]=u[1,1,2]-u[1,3,2]
      du[2,3,1]=u[2,1,2]-u[2,3,2]

    else         # Evaluate dv/dt
         
      q12_x = u[1,1,1]
      q12_y = u[2,1,1]
      q23_x = u[1,2,1]
      q23_y = u[2,2,1]
      q31_x = u[1,3,1]
      q31_y = u[2,3,1]

      dot12=(q12_x*q12_x+q12_y*q12_y)^-1.5
      dot23=(q23_x*q23_x+q23_y*q23_y)^-1.5
      dot31=(q31_x*q31_x+q31_y*q31_y)^-1.5
      
      aux12_x = q12_x*dot12
      aux12_y = q12_y*dot12
      aux23_x = q23_x*dot23
      aux23_y = q23_y*dot23
      aux31_x = q31_x*dot31
      aux31_y = q31_y*dot31

      dv1_x = Gm[2]*aux12_x - Gm[3]*aux31_x
      dv1_y = Gm[2]*aux12_y - Gm[3]*aux31_y
      dv2_x = Gm[3]*aux23_x - Gm[1]*aux12_x 
      dv2_y = Gm[3]*aux23_y - Gm[1]*aux12_y 
      dv3_x = Gm[1]*aux31_x - Gm[2]*aux23_x 
      dv3_y = Gm[1]*aux31_y - Gm[2]*aux23_y 
      
      du[1,1,2] = dv1_x
      du[2,1,2] = dv1_y
      du[1,2,2] = dv2_x
      du[2,2,2] = dv2_y
      du[1,3,2] = dv3_x 
      du[2,3,2] = dv3_y
      
    end # if

    return nothing

end






