function [ Ty,t1,t2,t3,t4 ] = Calculate_Ty_MN_at_Dx( i,j,Ex,Ey,Dx,Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy )

% G_y is at Hz
Ty=zeros(size(Ex));

c=299792458;
mu_o=4*pi*10^-7;

eps_o=(1/(c*c*mu_o));

% Ty=\px(Tyx) +\py(Tyy)
% Tyx= DyEx
% Tyy=-1/2DyEy

%% PLACE Ty at Ex (i+1/2, j), n


	Dy2=(0.5)*(Dy(i+1,j)+Dy(i+1,j-1));  % Dy moved to i+1,j-1/2, relative to Hz
	Ex2=(0.5)*(Ex(i,j)+Ex(i+1,j));		% Ex moved to i+1/2,j

	Dy1=(0.5)*(Dy(i,j)+Dy(i,j-1));		% Dy moved to i,j-1/2
	Ex1=(0.5)*(Ex(i,j)+Ex(i-1,j));		% Ex moved to i-1/2,j

% t1 \px (DyEx)

	t1=(1/dx)*(Dy2.*Ex2-Dy1.*Ex1);  % x derivative taken at i,j-1/2

% t2=(0.5)*\py (DxEx)

	t2=(0.5)*(1/(2*dy)*(Dx(i,j+1).*Ex(i,j+1)-Dx(i,j-1).*Ex(i,j-1)));

% t3=(0.5)*\py(DyEy)

% Dy2=(0.5)*(Dy(i,j)+Dy(i+1,j));
% Ey2=(0.5)*(Ey(i,j)+Ey(i+1,j));
 
% Dy1=(0.5)*(Dy(i,j-1)+Dy(i+1,j-1));
% Ey1=(0.5)*(Ey(i,j-1)+Ey(i+1,j-1));
 
% t3=(0.5)*(1/dy)*(Dy2.*Ey2-Dy1.*Ey1);

% Draw out yee grid to see why it has to be -2\
% The Ey component only meets a centered difference then .
	Dy2=(0.5)*(Dy(i,j+1)+Dy(i+1,j+1));
	Ey2=(0.5)*(Ey(i,j+1)+Ey(i+1,j+1));

	Dy1=(0.5)*(Dy(i,j-1)+Dy(i+1,j-2));
	Ey1=(0.5)*(Ey(i,j-2)+Ey(i+1,j-2));

	t3=(0.5)*(1/(3*dy))*(Dy2.*Ey2-Dy1.*Ey1);
	
% Apparently this was used in the working version ? why -2 on one side i dont know
	
	
	% Dy2=(0.5)*(Dy(i,j+1)+Dy(i+1,j+1));
	% Ey2=(0.5)*(Ey(i,j+1)+Ey(i+1,j+1));

	% Dy1=(0.5)*(Dy(i,j-2)+Dy(i+1,j-2));
	% Ey1=(0.5)*(Ey(i,j-2)+Ey(i+1,j-2));



% t5= \py(BzHz)

t4=(0.5)*(1/(3*dy))*(Hz(i,j+1).*Bz(i,j+1)-Hz(i,j-2).*Bz(i,j-2));

Ty(i,j)=-t1+t2-t3+t4;

end

