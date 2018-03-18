function [ Tx,t1,t2,t3,t4,t5 ] = calculate_Tx_AB_TE(i,j, Bx,By,Hx,Bx_n_prev,By_n_prev,Hx_n_prev,Hy_n_prev,Hy,Ez,Dz,dx,dy )


Tx=zeros(size(Bx));

% Tx= \px (Txx) +\py (Txy) +\pz (Txz);

%% Calculate at n, bring H to E temporally

Hy(i,j)=0.5.*(Hy(i,j)+Hy_n_prev(i,j));
By(i,j)=0.5*(By(i,j)+By_n_prev(i,j));

Hx(i,j)=0.5.*(Hx(i,j)+Hx_n_prev(i,j));
Bx(i,j)=0.5*(Bx(i,j)+Bx_n_prev(i,j));

% Tx=.5*(t1-t2+t3)-t4

% t1=\px(DzEz)

t1=.5*(1/dx).*(Dz(i,j).*Ez(i,j)-Dz(i-1,j).*Ez(i-1,j));

% t2=\px(BxHx)

Bx_avg(i,j)=0.5*(Bx(i,j)+Bx(i,j+1));
Hx_avg(i,j)=0.5*(Hx(i,j)+Hx(i,j+1));

t2=-.5*(1/dx).*(Bx_avg(i,j).*Hx_avg(i,j)-Bx_avg(i-1,j).*Hx_avg(i-1,j));

% t3=\px(ByHy)

t3=.5*(1/(2*dx)).*(By(i+1,j).*Hy(i+1,j)-By(i-1,j).*Hy(i-1,j));

% t4= \py(BxHy)
Bx_avg(i,j)=0.5*(Bx(i,j+1)+Bx(i-1,j+1));
Hy_avg(i,j)=0.5.*(Hy(i,j)+Hy(i,j+1));

t4=-.5*1*(1/dy).*( Bx_avg(i,j).*Hy_avg(i,j)-Bx_avg(i,j-1).*Hy_avg(i,j-1));

% t5 =\py(HxBy)

Hx_avg(i,j)=0.5*(Hx(i,j+1)+Hx(i-1,j+1));
By_avg(i,j)=0.5.*(By(i,j)+By(i,j+1));

t5=-.5*(1/dy).*( Hx_avg(i,j).*By_avg(i,j)-Hx_avg(i,j-1).*By_avg(i,j-1));


%Tx(i,j)=0.5.*((t1-t2+t3)-t4-t5);
Tx(i,j)=t1+t2+t3+t4+t5;

end

