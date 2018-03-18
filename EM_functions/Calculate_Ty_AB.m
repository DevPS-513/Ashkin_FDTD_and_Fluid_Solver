function [ Ty,t1,t2,t3,t4 ] = Calculate_Ty_AB( i,j,Ex,Ey,Dx,Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy )

Ty=zeros(size(Ex));


c=299792458;
mu_o=4*pi*10^-7;

eps_o=(1/(c*c*mu_o));


%
%% PLACE TY at Ex i+1/2, j
% t1=\px (EyEx)

Dy2=(0.5)*(Dy(i+1,j)+Dy(i+1,j-1));
Ex2=(0.5)*(Ex(i,j)+Ex(i+1,j));

Dy1=(0.5)*(Dy(i,j)+Dy(i,j-1));
Ex1=(0.5)*(Ex(i,j)+Ex(i-1,j));

t1=(1/dx)*(Dy2.*Ex2-Dy1.*Ex1);

% t2=eps_o*.5*\py(ExEx)


Ey2=(0.5)*(Ey(i+1,j)+Ey(i+1,j-1));
Dx2=(0.5)*(Dx(i,j)+Dx(i+1,j));

Ey1=(0.5)*(Ey(i,j)+Ey(i,j-1));
Dx1=(0.5)*(Dx(i,j)+Dx(i-1,j));

t2=(1/dx)*(Ey2.*Dx2-Ey1.*Dx1);

% t3=\px (DxEx)


t3=(1/(2*dy)*(Dx(i,j+1).*Ex(i,j+1)-Dx(i,j-1).*Ex(i,j-1)));

%t4=\py(DyDy)

Dy2=(0.5)*(Dy(i,j)+Dy(i+1,j));
Ey2=(0.5)*(Ey(i,j)+Ey(i+1,j));

Dy1=(0.5)*(Dy(i,j-1)+Dy(i+1,j-1));
Ey1=(0.5)*(Ey(i,j-1)+Ey(i+1,j-1));

t4=(1/dy)*(Ey2.*Dy2-Ey1.*Dy1);

% t5= \py(BzHz)

t5=(1/dy)*(Hz(i,j).*Bz(i,j)-Hz(i,j-1).*Bz(i,j-1));

Ty(i,j)=1/2*(-t1-t2+t3-t4+t5);

end

