function [ Ty,t1,t2,t3,t4 ] = Calculate_Ty_MN_v2( i,j,Ex,Ex_n_prev,Ey,Ey_n_prev,...
    Dx,Dx_n_prev,Dy,Dy_n_prev,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy )

% G_y is at Hz
Ty=zeros(size(Ex));

c=299792458;
mu_o=4*pi*10^-7;

eps_o=(1/(c*c*mu_o));

% Ty=\px(Tyx) +\py(Tyy)
% Tyx= DyEx
% Tyy=-1/2DyEy

%% PLACE Ty at Ex (i+1/2, j), n

% place Ty ay Ex in Space (i+1/2,j)
% Place Ty at Bz in time


Dy_avg=zeros(size(Ex));
Ex_avg=zeros(size(Ex));


Dy_avg=zeros(size(Ex));
Ex_avg=zeros(size(Ex));

Dy_avg(i,j)=.25*(Dy(i,j)+Dy(i,j-1)+Dy_n_prev(i,j)+Dy_n_prev(i,j-1));
Ex2=.25*(Ex(i,j)+Ex_n_prev(i,j)+Ex(i+1,j)+Ex_n_prev(i+1,j));
Ex1=.25*(Ex(i,j)+Ex_n_prev(i,j)+Ex(i-1,j)+Ex_n_prev(i-1,j));

	
t1=-1*(1/dx)*(Dy_avg(i+1,j).*Ex2-Dy_avg(i,j).*Ex1);

% t1 \px (DyEx)
% t2=(0.5)*\py (DxEx) x derivative about i+12, j
Dx_at_Bt=zeros(size(Ex));
Ex_at_Bt=zeros(size(Ex));


Dx_at_Bt(i,j)=0.5*(Dx(i,j)+Dx_n_prev(i,j));
Ex_at_Bt(i,j)=0.5*(Ex(i,j)+Ex_n_prev(i,j));

t2=(0.5)*(1/(2*dy)*(Dx_at_Bt(i,j+1).*Ex_at_Bt(i,j+1)-Dx_at_Bt(i,j-1).*Ex_at_Bt(i,j-1)));

% t3=(0.5)*\py(DyEy)

% Dy2=(0.5)*(Dy(i,j)+Dy(i+1,j));
% Ey2=(0.5)*(Ey(i,j)+Ey(i+1,j));
 
% Dy1=(0.5)*(Dy(i,j-1)+Dywwww(i+1,j-1));
% Ey1=(0.5)*(Ey(i,j-1)+Ey(i+1,j-1));
 
% t3=(0.5)*(1/dy)*(Dy2.*Ey2-Dy1.*Ey1);

% Draw out yee grid to see why it has to be -2\
% The Ey component only meets a centered difference then .

Dy2=zeros(size(Ex));
Ey2=zeros(size(Ex));

Dy1=zeros(size(Ex));
Ey1=zeros(size(Ex));

	Dy2(i,j)=(0.25)*(Dy(i,j)+Dy(i+1,j)+Dy_n_prev(i,j)+Dy_n_prev(i+1,j));
	Ey2(i,j)=(0.25)*(Ey(i,j)+Ey(i+1,j)+Ey_n_prev(i,j)+Ey_n_prev(i+1,j));

Dy1(i,j)=(0.25)*(Dy(i,j-1)+Dy(i+1,j-1)+Dy_n_prev(i,j-1)+Dy_n_prev(i+1,j-1));
Ey1(i,j)=(0.25)*(Ey(i,j-1)+Ey(i+1,j-1)+Ey_n_prev(i,j-1)+Ey_n_prev(i+1,j-1));


	t3=-(0.5)*(1/(dy))*(Dy2(i,j).*Ey2(i,j)-Dy1(i,j).*Ey1(i,j));
	
% Apparently this was used in the working version ? why -2 on one side i dont know
	
	
	% Dy2=(0.5)*(Dy(i,j+1)+Dy(i+1,j+1));
	% Ey2=(0.5)*(Ey(i,j+1)+Ey(i+1,j+1));

	% Dy1=(0.5)*(Dy(i,j-2)+Dy(i+1,j-2));
	% Ey1=(0.5)*(Ey(i,j-2)+Ey(i+1,j-2));



% t5= \py(BzHz)

t4=(0.5)*(1/(dy))*(Hz(i,j).*Bz(i,j)-Hz(i,j-1).*Bz(i,j-1));

Ty(i,j)=t1+t2+t3+t4;



end

