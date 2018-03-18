function [ Dx,Dy ] = diffuse_visc( u,v,dx,dy,mu )
%% i


%[Nx,kNy]=size(u);

[Nxu, Nyu]=size(u);
[Nxv,Nyv]=size(v);
Dx=zeros(Nxu,Nyu);
Dy=zeros(Nxv,Nyv);

%% diffuse  the x-component

% index(i,j) refers to u(i+1/2,j) and v(i,j+1/2)

for q=2:Nxu-1
    
    for k=2:Nyu-1
        
   
        mu_av= 0.25*(mu(q+1,k)+mu(q+1,k+1)+mu(q,k+1)+mu(q,k)); %mu (i+.5,j+.5)
        
         mu_ym= 0.25*(mu(q,k)+mu(q+1,k)+mu(q,k-1)+mu(q+1,k-1)); %mu (i+.5,j-.5)
        

        T_1=1/dx*( 2*mu(q+1,k)*(u(q+1,k)-u(q,k))/dx-2*mu(q,k)*(u(q,k)-u(q-1,k))/dx);
        T_2=1/dy*(mu_av*( (u(q,k+1)-u(q,k))/dy+(v(q+1,k)-v(q,k))/dx));
        T_3=-1*mu_ym*((u(q,k)-u(q,k-1))/dy +(v(q+1,k-1)-v(q,k-1))/dx);
        
        Dx(q,k)=T_1+T_2+T_3;



    end
end
   
    %% diffuse Ay component
 for q=2:Nxv-1
    
    for k=2:Nyv-1

           mu_av= 0.25*(mu(q+1,k)+mu(q+1,k+1)+mu(q,k+1)+mu(q,k)); %mu (i+.5,j+.5)
        
         mu_xm= 0.25*(mu(q,k)+mu(q,k+1)+mu(q-1,k+1)+mu(q-1,k)); %mu (i+.5,j-.5)
        
        
        
 T_1=(1/dx)*(mu_av*( (u(q,k+1)-u(q,k))/dy +(v(q+1,k)-v(q,k))/dx));
 T_2= -(1/dx)*mu_xm*((u(q-1,k+1)-u(q-1,k))/dy+(v(q,k)-v(q-1,k))/dx);
 T_3=(1/dy)*(2*mu(q,k+1)*(v(q,k+1)-v(q,k))/dy-2*mu(q,k)*(v(q,k)-v(q,k-1))/dy);

 Dy(q,k)=T_1+T_2+T_3;

    end
    
end

end

