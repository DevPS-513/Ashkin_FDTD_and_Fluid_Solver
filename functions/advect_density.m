function [ rho ] = advect_density( u,v,dx,dy,dt,rho,mu_o )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



[Nxr,Nyr]=size(rho);

r=rho;

for q=2:Nxr-1
    for k=2:Nyr-1
        
        
        r1=0.5*(r(q+1,k)+r(q,k));
        r2=0.5*(r(q-1,k)+r(q,k));
        r3=0.5*(r(q,k+1)+r(q,k));
        r4=0.5*(r(q,k-1)+r(q,k));
        
        
        Div_x=(1/dx)*(r1*u(q,k)-r2*u(q-1,k));
        Div_y=(1/dy)*(r3*v(q,k)-r4*v(q,k-1));
        Div=Div_x+Div_y;
        
        Lap=(1/(dx^2))*(r(q+1,k)-2*r(q,k)+r(q-1,k))+(1/(dy^2))*(r(q,k+1)-2*r(q,k)+r(q,k-1));
       
        rho(q,k)=r(q,k)-dt*(Div)+dt*mu_o*Lap;
        
        
    end
    
    
end


end

