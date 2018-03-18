function [ Ty,t1,t2,t3,t4,t5 ] = calculate_Ty_AB_TE(i,j, Bx,By,Hx,Bx_n_prev,By_n_prev,Hx_n_prev,Hy_n_prev,Hy,Ez,Dz,dx,dy )



Ty=zeros(size(Ez));
% Ty=\px(TyTx) + \py( Tyy) +\pz(Tzy)
By_avg=zeros(size(By));
Bx_avg=zeros(size(By));
Hx_avg=zeros(size(Hx));
Hy_avg=zeros(size(Hy));
% t1=\px(ByHx)
% we place Ty at Hx (i+1/2,j)


By_avg(i,j)=.25*(By(i,j)+By(i,j-1)+By_n_prev(i,j)+By_n_prev(i,j-1));
Hx_avg(i,j)=.25*(Hx(i,j)+Hx_n_prev(i,j)+Hx(i-1,j)+Hx_n_prev(i-1,j));


t1=-.5*(1/dx)*(By_avg(i+1,j).*Hx_avg(i+1,j)-By_avg(i,j).*Hx_avg(i,j));

%t2=\px HyBx

Hy_avg(i,j)=.25*(Hy(i,j)+Hy(i,j-1)+Hy_n_prev(i,j)+Hy_n_prev(i,j-1));
Bx_avg(i,j)=.25*(Bx(i,j)+Bx_n_prev(i,j)+Bx(i-1,j)+Bx_n_prev(i-1,j));


t2=-.5*(1/dx)*(Hy_avg(i+1,j).*Bx_avg(i+1,j)-Hy_avg(i,j).*Bx_avg(i,j));

% t3=\py DzEz

t3=.5*(1/dy)*(Dz(i,j).*Ez(i,j)-Dz(i,j-1).*Ez(i,j-1));

% t4= \py BxHx

t4=.5*(1/(2*dy))*(Bx(i,j+1).*Hx(i,j+1)-Bx(i,j-1).*Hx(i,j-1));

% t= \py(ByHy)

By_avg(i,j)=.5*(By(i,j)+By(i+1,j));
Hy_avg(i,j)=.5*(Hy(i,j)+Hy(i+1,j));

t5=-1*(1/dy)*(By_avg(i,j).*Hy_avg(i,j)-By_avg(i,j-1).*Hy_avg(i,j-1));




%Ty(i,j)=0.5*(-t1-t2+t3+t4-t5);
Ty(i,j)=t1+t2+t3+t4+t5;






end

