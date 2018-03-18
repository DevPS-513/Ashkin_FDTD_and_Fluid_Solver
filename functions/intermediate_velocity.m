function [ u_star,v_star ] = intermediate_velocity( u,v,rho,gx,gy,fx,fy,Ax,Ay,Dx,Dy,dt )
%UNTITLED5 Summary of this function goes here
%   Detailed ekxplanation goes here

[Nxu, Nyu]=size(u);
[Nxv,Nyv]=size(v);
u_star=zeros(Nxu,Nyu);
v_star=zeros(Nxv,Nyv);

for q=2:Nxu-1
    for k=2:Nyu-1
        
        rhox=1/2*(rho(q+1,k)+rho(q,k));
        
        u_star(q,k)=u(q,k)+dt*(-1*Ax(q,k)+gx(q,k)+(1/rhox)*(Dx(q,k)+fx(q,k)));
        
        
    end
    
end
    
    
 for q=2:Nxv-1
    for k=2:Nyv-1
        
        rhoy=1/2*(rho(q,k+1)+rho(q,k));
        v_star(q,k)=v(q,k)+dt*(-1*Ay(q,k)+gy(q,k)+(1/rhoy)*(Dy(q,k)+fy(q,k)));
        
        
    end
    
end


end

