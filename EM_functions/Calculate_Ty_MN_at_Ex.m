function [ Ty ] = Calculate_Ty_MN_at_Ex( i,j,Ex,Ey,Dx,Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy )

Ty=zeros(size(Ex));



% div(T)=T_1+T_2+T_3
% i,j will refer to Ex location


t1=(1/dx)*(Dy(i+1,j).*Ex2-Dy(i,j).*Ex1);

% t2= .5(1/dy) DxEx

t2=0.5*(1/dy)*(Dx(i,j+1).*Ex(i,j+1)-Dx(i,j).*Ex(i,j));

% t3= .5*\py(DyEy)
Dy2=(0.25)*(Dy(i,j)+Dy(i+1,j)+Dy(i,j+1)+Dy(i+1,j+1));
Dy1=(0.25)*(Dy(i,j)+Dy(i+1,j)+Dy(i,j-1)+Dy(i+1,j-1));

Ey2=(0.25)*(Ey(i,j)+Ey(i+1,j)+Ey(i,j+1)+Ey(i+1,j+1));
Ey1=(0.25)*(Ey(i,j)+Ey(i+1,j)+Ey(i,j-1)+Ey(i+1,j-1));

t3=(0.5*(1/dy)*(Dy2.*Ey2-Dy1.*Ey1));

% t4= (0.5)*\py(BzHz)

t4=(0.5)*(1/(2*dy))*(Bz(i,j+1).*Hz(i,j+1)-Bz(i,j-1).*Hz(i,j-1));

% t5=\py(DyEy)
t5=2*t3;


Ty=-t1+t2+t3+t4-t5;


end

