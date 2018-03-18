function [ Tx,t1,t2,t3,t4,t5 ] = calculate_Tx_MN_TE(i,j, Bx,By,Hx,Bx_n_prev,By_n_prev,Hx_n_prev,Hy_n_prev,Hy,Ez,Dz,dx,dy )


Tx=zeros(size(Bx));


%% Calculate at n, bring H to E temporally


% Tx=.5*(t1-t2+t3)-t4

% t1=\px(DzEz)

t1=.5*(1/dx).*(Dz(i,j).*Ez(i,j)-Dz(i-1,j).*Ez(i-1,j));

% t2=\px(BxHx)

Bx_avg(i,j)=0.5*(Bx(i,j)+Bx(i,j+1));
Hx_avg(i,j)=0.5*(Hx(i,j)+Hx(i,j+1));

t2=-.5*(1/dx).*(Bx_avg(i,j).*Hx_avg(i,j)-Bx_avg(i-1,j).*Hx_avg(i-1,j));

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
Tx(i,j)=t1+t2+t3+t4;


end

