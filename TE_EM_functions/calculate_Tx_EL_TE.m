function [ Tx,t1,t2,t3,t4,t5 ] = calculate_Tx_EL_TE(i,j, Bx,By,Hx,Bx_n_prev,By_n_prev,Hx_n_prev,Hy_n_prev,Hy,Ez,Dz,dx,dy )

meters=1;
nm=meters*1e-9;
fs=1e-15;
mu_o=4*pi*10^-7;
c=299792458;
eps_o=(1/(c*c*mu_o));
eta_o=sqrt(mu_o/eps_o);
Kg=1;

Tx=zeros(size(Bx));


t1=.5*eps_o.*(1/dx).*(Ez(i,j).*Ez(i,j)-Ez(i-1,j).*Ez(i-1,j));

% t2=\px(BxHx)

Bx_avg(i,j)=0.5*(Bx(i,j)+Bx(i,j+1));

t2=-.5*(1/dx).*(Bx_avg(i,j).*Bx_avg(i,j)-Bx_avg(i-1,j).*Bx_avg(i-1,j))./mu_o;

% t3=\px(ByHy)

Hy=0.5.*(Hy+Hy_n_prev);
By=0.5*(By+By_n_prev);


t3=.5*(1/(2*dx)).*(By(i+1,j).*Hy(i+1,j)-By(i-1,j).*Hy(i-1,j));

% t4= \py(BxHy)
Bx_avg(i,j)=0.5*(Bx(i,j+1)+Bx(i-1,j+1));
Hy_avg(i,j)=0.5.*(Hy(i,j)+Hy(i,j+1));

t4=-1*(1/dy).*( Bx_avg(i,j).*Hy_avg(i,j)-Bx_avg(i,j-1).*Hy_avg(i,j-1));

t5=0;

%Tx(i,j)=0.5.*(t1-t2+t3)-t4;
Tx(i,j)=t1+t2+t3+t4+t5;

end

