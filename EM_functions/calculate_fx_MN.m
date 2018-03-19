function [ fx, fy] = calculate_fx_MN(i,j,Hz,Hz_n_prev,Bz,Bz_n_prev,Ex,Ey,Dx,Dy,dx,dy )




% f_MN=-1/2[-EX(nabXD)+DX(\nabXE)-HX(\nabXB)+BX(\nabXH)];
% f_MN=-1/2*(-f1m+f2m-f3m+f4m)

% f1m= EX(\na bX D)
% Place fx at Ey
Dx_diff_y2=(Dx(i,j+1)-Dx(i,j));
Dx_diff_y1=(Dx(i-1,j+1)-Dx(i-1,j));
Dx_diff=1/2*(Dx_diff_y2+Dx_diff_y1);

f1x=Ey(i,j).*((1/(2*dx))*(Dy(i+1,j)-Dy(i-1,j))-(1/(dy)).*(Dx_diff));

% f2m= DX(\nab X E)
% by symmetery sceond term is simply first term but swap E and D
Ex_diff_y2=(Ex(i,j+1)-Ex(i,j));
Ex_diff_y1=(Ex(i-1,j+1)-Ex(i-1,j));
Ex_diff=1/2*(Ex_diff_y2+Ex_diff_y1);


f2x=Dy(i,j).*((1/(2*dx))*(Ey(i+1,j)-Ey(i-1,j))-(1/(dy)).*(Dx_diff));
% f3m=HX(nab X B )
Hz_avg=(1/4)*(Hz(i,j)+Hz(i-1,j)+Hz_n_prev(i,j)+Hz_n_prev(i-1,j));
Bz_x2=(1/2).*(Bz(i,j)+Bz_n_prev(i,j));
Bz_x1=(1/2).*(Bz(i-1,j)+Bz_n_prev(i-1,j));

f3x=1*Hz_avg.*(1/dx).*(Bz_x2-Bz_x1);


% f4m= BX(nabXH)
Bz_avg=(1/4)*(Bz(i,j)+Bz(i-1,j)+Bz_n_prev(i,j)+Bz_n_prev(i-1,j));
Hz_x2=(1/2).*(Hz(i,j)+Hz_n_prev(i,j));
Hz_x1=(1/2).*(Hz(i-1,j)+Hz_n_prev(i-1,j));

f4x=Bz_avg.*(1/dx).*(Hz_x2-Hz_x1);

fx=1*(1/2).*(-1.*f1x+f2x-1.*f3x+f4x);
fy=0;
end

