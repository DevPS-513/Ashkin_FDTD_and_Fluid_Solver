function [ Tx ,t1,t2,t3,t4,t5,t6] = Calculate_Tx_AB_at_Ey( i,j,Ex,Ey,Dx,Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy )

Tx=zeros(size(Ex));
% bring Ey to Ex locations
% being Dy to Dx locations

ic=[i(1)-2:i(end)+2];
jc=[j(1)-2:j(end)+2];

% Initialize avg matricies


Ex_av=zeros(size(Ex));

Dx_av=zeros(size(Ex));

Hz_av=zeros(size(Hz));
Bz_av=zeros(size(Bz));

% creat avg matricies

% T1=t1+t2+t3-t4




Ex_av(ic,jc)=(1/2)*( Ex(ic,jc)+Ex(ic,jc+1) );
Dx_av(ic,jc)=(1/2)*(Dx(ic,jc)+Dx(ic,jc+1));

Hz_av(ic,jc)=(1/2)*(Hz(ic,jc)+Hz_n_prev(ic,jc));
Bz_av(ic,jc)=(1/2)*(Bz(ic,jc)+Bz_n_prev(ic,jc));


% t1=\px(DxEx)

t1=(1/(3*dx)).*(Dx_av(i+1,j).*Ex_av(i+1,j)-Dx_av(i-2,j).*Ex(i-2,j));

% t2=\px DyEy

t2=(1/(2*dx)).*(Dy(i+1,j).*Ey(i+1,j)-Dy(i-1,j).*Ey(i-1,j));

% t3=\px(BzHz) ( double diff for now? )
t3=(1/(3*dx)).*(Bz_av(i+1,j).*Hz_av(i+1,j)-Bz_av(i-2,j).*Hz_av(i-2,j));

% t4=2*DxEx
t4=(2/(3*dx)).*(Dx_av(i+1,j).*Ex_av(i+1,j)-Dx_av(i-2,j).*Ex(i-2,j));

T1=t1+t2+t3-t4;

% requires more averageing
Dx_2=(1/2).*(Dx(i,j+1)+Dx(i-1,j+1));
Ey_2=(1/2).*(Ey(i,j)+Ey(i,j+1));

Dx_1=(1/2).*(Dx(i,j)+Dx(i-1,j));
Ey_1=(1/2).*(Ey(i,j)+Ey(i,j-1));

t5=(1/dy).*(Dx_2.*Ey_2  -  Dx_1.*Ey_1);

% requires more averagine
Ex_2=(1/2).*(Ex(i,j+1)+Ex(i-1,j+1));
Dy_2=(1/2).*(Dy(i,j)+Dy(i,j+1));

Ex_1=(1/2).*(Ex(i,j)+Ex(i-1,j));
Dy_1=(1/2).*(Dy(i,j)+Dy(i,j-1));

t6=(1/dy)*(Ex_2.*Dy_2-Ex_1.*Dy_1);

T2=-t5-t6;

Tx(i,j)=1/2*(T1+T2);
end

