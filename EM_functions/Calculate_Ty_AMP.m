function [ Ty,t1,t2,t3,t4 ] = Calculate_Ty_AMP( i,j,Ex,Ey,Dx,Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy )

Ty=zeros(size(Ex));


c=299792458;
mu_o=4*pi*10^-7;

eps_o=(1/(c*c*mu_o));


%% PLACE TY at Ex i+1/2, j
% t1=\px (EyEx)

Ey2=(0.5)*(Ey(i+1,j)+Ey(i+1,j-1));
Ex2=(0.5)*(Ex(i,j)+Ex(i+1,j));

Ey1=(0.5)*(Ey(i,j)+Ey(i,j-1));
Ex1=(0.5)*(Ex(i,j)+Ex(i-1,j));

t1=eps_o*(1/dx)*(Ey2.*Ex2-Ey1.*Ex1);

% t2=eps_o*.5*\py(ExEx)


t2=(0.5)*eps_o*(1/(2*dy)*(Ex(i,j+1).^2-Ex(i,j-1).^2));


%t3=eps_o*.5*\py(EyEy)

Ey2=(0.5)*(Ey(i,j)+Ey(i+1,j));
Ey1=(0.5)*(Ey(i,j-1)+Ey(i+1,j-1));

t3=(0.5)*eps_o*(1/dy)*(Ey2.^2-Ey1.^2);


% t4=(0.5)*mu_o*\py(HzHz)

t4=(0.5)*(1/mu_o)*(1/dy)*(Bz(i,j).^2-Bz(i,j-1).^2);

Ty(i,j)=-t1+t2-t3+t4;


end

