function [ u,v ] = velocity_correct(u,v,u_star,v_star,dx,dy,dt,rho,p)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goeks here


[Nxu, Nyu]=size(u);
[Nxv,Nyv]=size(v);


for q=2:Nxu-1
    
    for k=2:Nyu-1
        
        C1=dt/(dx*.5*(rho(q+1,k)+rho(q,k)));
        
        u(q,k)=u_star(q,k)-C1*(p(q+1,k)-p(q,k));
        
        
        
        
        
        
    end
    
    
end



for q=2:Nxv-1
    
    for k=2:Nyv-1
        
        C2=dt/(dy*.5*(rho(q,k+1)+rho(q,k)));
        
        v(q,k)=v_star(q,k)-C2*(p(q,k+1)-p(q,k));
        
        
        
        
        
        
    end
    
    
end



end

