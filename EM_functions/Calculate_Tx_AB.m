function [ Tx,Ey_av,t1,t2,t3,t4,T2,T3 ] = Calculate_Tx_AB( i,j,Ex,Ey,Dx,Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy )

Tx=zeros(size(Ex));
% bring Ey to Ex locations
% being Dy to Dx locations

ic=i;
jc=j;

% Initialize avg matricies


Ey_av=zeros(size(Ex));
Dy_av=zeros(size(Ex));

Hz_av=zeros(size(Ex));
Bz_av=zeros(size(Ex));

% creat avg matricies



Ey_av(ic,jc)=1/4*(Ey(ic,jc)+Ey(ic+1,jc)+Ey(ic,jc-1)+Ey(ic+1,jc-1));
Dy_av(ic,jc)=1/4*(Dy(ic,jc)+Dy(ic+1,jc)+Dy(ic,jc-1)+Dy(ic+1,jc-1));

Hz_av(ic,jc)=1/4*(Hz(ic,jc)+ Hz(ic,jc-1)+ Hz_n_prev(ic,jc)+ Hz_n_prev(ic,jc-1));
Bz_av(ic,jc)=1/4*(Bz(ic,jc)+ Bz(ic,jc-1)+ Bz_n_prev(ic,jc)+ Bz_n_prev(ic,jc-1));


% div(T)=T_1+T_2+T_3
% i,j will refer to Ex location

%% T1= \px(Txx)=t1+t2+t3-t4

t1=(1/2)*(1/(2*dx)).*(Dx(i+1,j).*Ex(i+1,j)-Dx(i-1,j).*Ex(i-1,j));

% t2=(1/2)\px(DxEy)

t2=(1/2)*(1/(2*dx)).*(Dy_av(i+1,j).*Ey_av(i+1,j)-Dy_av(i-1,j).*Ey_av(i-1,j));

% t3=(1/2)\px(BzHz)

t3=(1/2)*(1/(2*dx)).*(Bz_av(i+1,j).*Hz_av(i+1,j)-Bz_av(i-1,j).*Hz_av(i-1,j));

% t4=\px(DxEx)

t4=(1/(2*dx)).*(Dx(i+1,j).*Ex(i+1,j)-Dx(i-1,j).*Ex(i-1,j));

T1=t1+t2+t3-t4;

%% T2=-DxEy

T2=-1*(1/(2*dy)).*(Dx(i,j+1).*Ey_av(i,j+1)-Dx(i,j-1).*Ey_av(i,j-1));

T3=-1*(1/(2*dy)).*(Ex(i,j+1).*Dy_av(i,j+1)-Ex(i,j-1).*Dy_av(i,j-1));


Tx(i,j)=T1+1/2*(T2+T3);

%% T


end

