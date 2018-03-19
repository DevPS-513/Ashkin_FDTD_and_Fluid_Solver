function [ fy, fJy,fpy ] = calculate_fy_MS( i,j,Hz,Hz_n_prev,Jmz,Jmz_n_prev,Ex,Ey,Jex,Jey,Px,Py,Mzs,dx,dy )

mu_o=4*pi*10^-7;
c=299792458;
eps_o=(1/(c*c*mu_o));
% i=[NPML_x:1:Nx-NPML_x];
% j=[NPML_y:1:Ny-NPML_y];


% f= (P*nab)E+(M*nab)H+(\ptP)X(mu_oH)-(\ptM)X(\eps_o E)


% place fy ontop of Ey and Px i j will refer to i+1/2,j

%% P TERM
Px_av=1/4*(Px(i,j)+Px(i-1,j)+Px(i,j+1)+Px(i-1,j+1));

fpy=Px_av*(1/(2*dx)).*(Ey(i+1,j)-Ey(i-1,j))+Py(i,j)*(1/(2*dy)).*(Ey(i,j+1)-Ey(i,j-1));

%% JX term
Jex_av=1/4*(Jex(i,j)+Jex(i-1,j)+Jex(i,j+1)+Jex(i-1,j+1));
Hz_av=1/4*( Hz(i,j)+ Hz(i-1,j) +Hz_n_prev(i,j)+ Hz_n_prev(i-1,j) );

fJy=Jex_av.*mu_o.*Hz_av;
%% JM term
Jmz_av=1/4*( Jmz(i,j)+ Jmz(i-1,j) +Jmz_n_prev(i,j)+ Jmz_n_prev(i-1,j) );
Ex_av=1/4*(Ex(i,j)+Ex(i-1,j)+Ex(i,j+1)+Ex(i-1,j+1));

fJmy=Jmz_av.*eps_o.*Ex_av;


% other terms are zero for a TM mode.
 fy(i,j)=fpy+fJy+fJmy;

end

