function [ Ax,Ay ] = advect( u,v,dx,dy )


%[Nx,Ny]=sikze(u);

[Nxu, Nyu]=size(u);
[Nxv,Nyv]=size(v);
Ax=zeros(Nxu,Nyu);
Ay=zeros(Nxv,Nyv);

%% Advect the x-component

% index(i,j) refers to u(i+1/2,j) and v(i,j+1/2)

for q=2:Nxu-1
    
    for k=2:Nyu-1


A_x1=((u(q+1,k)+u(q,k))/2)^2-((u(q,k)+u(q-1,k))/2)^2;
A_x2=((u(q,k+1)+u(q,k))*(v(q+1,k)+v(q,k)))/4-(( u(q,k)+u(q,k-1) )*(v(q+1,k-1)+v(q,k-1)) )/4;

Ax(q,k)=(1/(dx*dy))*(dy*A_x1+dx*A_x2);

    end
    
end   
    %% Advect y component
for q=2:Nxv-1
    
    for k=2:Nyv-1



A_y2=((v(q,k+1)+v(q,k))/2)^2-((v(q,k)+v(q,k-1))/2)^2;




A_y1=((v(q+1,k)+v(q,k))*(u(q,k+1)+u(q,k)))/4-(( v(q,k)+v(q-1,k) )*(u(q-1,k+1)+u(q-1,k)) )/4;


Ay(q,k)=(1/dx)*A_y1+(1/dy)*A_y2;

    end
    
end

end

