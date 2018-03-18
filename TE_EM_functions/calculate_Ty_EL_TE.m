function [ Ty,t1,t2,t3,t4,t5 ] = calculate_Ty_EL_TE(i,j, Bx,By,Hx,Bx_n_prev,By_n_prev,Hx_n_prev,Hy_n_prev,Hy,Ez,Dz,dx,dy )

meters=1;
nm=meters*1e-9;
fs=1e-15;
mu_o=4*pi*10^-7;
c=299792458;
eps_o=(1/(c*c*mu_o));
eta_o=sqrt(mu_o/eps_o);
Kg=1;


Ty=zeros(size(Ez));

By_avg=zeros(size(By));
Hx_avg=zeros(size(Hx));
Hy_avg=zeros(size(Hy));
% t1=\px(ByHx)
% we place Ty at Hx (i+1/2,j)


By_avg(i,j)=.25*(By(i,j)+By(i,j-1)+By_n_prev(i,j)+By_n_prev(i,j-1));
Hx_avg(i,j)=.25*(Hx(i,j)+Hx_n_prev(i,j)+Hx(i-1,j)+Hx_n_prev(i-1,j));


t1=-1*(1/dx)*(By_avg(i+1,j).*Hx_avg(i+1,j)-By_avg(i,j).*Hx_avg(i,j));


% t2=\py DzEz

t2=.5*eps_o*(1/dy)*(Ez(i,j).*Ez(i,j)-Ez(i,j-1).*Ez(i,j-1));

% t3= \py BxHx

t3=.5*(1/(2*dy))*(Hx(i,j+1).*Hx(i,j+1)-Hx(i,j-1).*Hx(i,j-1))*mu_o;

% t= \py(ByHy)

By_avg(i,j)=.5*(By(i,j)+By(i+1,j));
Hy_avg(i,j)=.5*(Hy(i,j)+Hy(i+1,j));

t4=-.5*(1/dy)*(Hy_avg(i,j).*Hy_avg(i,j)-Hy_avg(i,j-1).*Hy_avg(i,j-1))*mu_o;

t5=0;


%Ty(i,j)=-t1+0.5*(t2+t3-t4);
Ty(i,j)=t1+t2+t3+t4+t5;






end

