function [ Dx,Dy ] = diffuse( u,v,dx,dy,mu_o )


%[Nx,kNy]=size(u);

[Nxu, Nyu]=size(u);
[Nxv,Nyv]=size(v);
Dx=zeros(Nxu,Nyu);
Dy=zeros(Nxv,Nyv);

%% diffuse  the x-component

% index(i,j) refers to u(i+1/2,j) and v(i,j+1/2)

for q=2:Nxu-1
    
    for k=2:Nyu-1

D_1=(u(q+1,k)-2*u(q,k)+u(q-1,k))/(dx^2);
D_2=(u(q,k+1)-2*u(q,k)+u(q,k-1))/(dy^2);

Dx(q,k)=mu_o*(D_1+D_2);


    end
end
   
    %% diffuse Ay component
 for q=2:Nxv-1
    
    for k=2:Nyv-1
D_1y=(v(q+1,k)-2*v(q,k)+v(q-1,k))/(dx^2);
D_2y=(v(q,k+1)-2*v(q,k)+v(q,k-1))/(dy^2);

Dy(q,k)=mu_o*(D_1y+D_2y);


    end
    
end

end

