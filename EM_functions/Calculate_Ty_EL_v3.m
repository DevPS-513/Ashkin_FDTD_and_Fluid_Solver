function [ Ty,t1,t2,t3,t4 ] = Calculate_Ty_EL_v3( i,j,Ex,Ex_n_prev,Ey,Ey_n_prev,Dx,Dx_n_prev,...
    Dy,Dy_n_prev,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy )

Ty=zeros(size(Ex));


c=299792458;
mu_o=4*pi*10^-7;

eps_o=(1/(c*c*mu_o));



% place Ty ay Ex in Space (i+1/2,j)
% Place Ty at Bz in time


Dy2=.5* ( .25*(Dy(i+1,j)+Dy(i+2,j)+Dy(i+1,j-1)+Dy(i+2,j-1))...
    +.25*(Dy_n_prev(i+1,j)+Dy_n_prev(i+2,j)+Dy_n_prev(i+1,j-1)+Dy_n_prev(i+2,j-1)) );


Dy1=.5* ( .25*(Dy(i,j)+Dy(i-1,j)+Dy(i-1,j-1)+Dy(i-1,j))...
    +.25*(Dy_n_prev(i,j)+Dy_n_prev(i-1,j)+Dy_n_prev(i-1,j-1)+Dy_n_prev(i-1,j)) );


Ex2=.25*(Ex(i,j)+Ex_n_prev(i,j)+Ex(i+1,j)+Ex_n_prev(i+1,j));
Ex1=.25*(Ex(i,j)+Ex_n_prev(i,j)+Ex(i-1,j)+Ex_n_prev(i-1,j));

% t1=\px(DyEx)
t1=-(1/dx)*(Dy2.*Ex2-Dy1.*Ex1);

%t2=(1/2)\py(ExEx)

Ex_at_Bt=zeros(size(Ex));
Ex_at_Bt(i,j)=0.5*(Ex(i,j)+Ex_n_prev(i,j));

t2=(0.5)*(1/(2*dy)*(Ex_at_Bt(i,j+1).*Ex_at_Bt(i,j+1)-Ex_at_Bt(i,j-1).*Ex_at_Bt(i,j-1)));


% t3=(1/2)*eps_o\py(EyEy)




	Ey2=(0.25)*(Ey(i,j)+Ey(i+1,j)+Ey_n_prev(i,j)+Ey_n_prev(i+1,j));

Ey1=(0.25)*(Ey(i,j-1)+Ey(i+1,j-1)+Ey_n_prev(i,j-1)+Ey_n_prev(i+1,j-1));



t3=(0.5)*eps_o*(1/dy)*(Ey2.^2-Ey1.^2);


% t4=(0.5)*mu_o*\py(HzHz)

t4=(0.5)*mu_o*(1/dy)*(Hz(i,j).^2-Hz(i,j-1).^2);

%t5=\py(DyEy)
Dy2=(0.5)*(Dy(i,j)+Dy(i+1,j));
Ey2=(0.5)*(Ey(i,j)+Ey(i+1,j));

Dy1=(0.5)*(Dy(i,j-1)+Dy(i+1,j-1));
Ey1=(0.5)*(Ey(i,j-1)+Ey(i+1,j-1));

t5=(1/dy)*(Ey2.*Dy2-Ey1.*Dy1);


Ty(i,j)=t1+t2+t3+t4-t5;

end

