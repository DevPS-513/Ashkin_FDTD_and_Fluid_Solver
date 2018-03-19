function [ fx, fJx, fpx  ] = calculate_fx_AB( i,j,Hz,Hz_n_prev,Jmz,Jmz_n_prev_fr,Ex,Ey,Jex,Jey,Px,Py,Mzs,dx,dy )

mu_o=4*pi*10^-7;
c=299792458;
eps_o=(1/(c*c*mu_o));
% i=[NPML_x:1:Nx-NPML_x];
% j=[NPML_y:1:Ny-NPML_y];


% f=
% current density term
Hz_fr_x=(1/4)*(Hz(i,j)+Hz(i-1,j)+Hz_n_prev(i,j)+Hz_n_prev(i-1,j));
Jmz_fr_x=(1/4)*(Jmz(i,j)+Jmz(i-1,j)+Jmz_n_prev_fr(i,j)+Jmz_n_prev_fr(i-1,j));


fJx(i,j)=mu_o*Jey(i,j).*Hz_fr_x+eps_o.*Jmz_fr_x.*Ey(i,j);

% fpx-x
  
Ex_av_1=1/2*(Ex(i-1,j)+Ex(i-1,j+1));
Ex_av_2=1/2*(Ex(i,j)+Ex(i,j+1));
Px_av=1/4*(Px(i,j)+Px(i-1,j)+Px(i-1,j-1)+Px(i,j+1));
            
fpx(i,j)=Px_av.*(1/dx).*(Ex_av_2-Ex_av_1);

% fpy-x
Ex_av_y2=1/2*(Ex(i,j+1)+Ex(i-1,j+1));
Ex_av_y1=1/2*(Ex(i,j)+Ex(i-1,j));

fpx(i,j)=fpx(i,j)+Py(i,j).*(1/dy).*(Ex_av_y2-Ex_av_y1 );

fx=fJx(i,j)+fpx(i,j);
  




end

