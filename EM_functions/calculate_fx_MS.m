function [ fx, fJex, fpx  ] = calculate_fx_MS( i,j,Hz,Hz_n_prev,Jmz,Jmz_n_prev,Ex,Ey,Jex,Jey,Px,Py,Mzs,dx,dy )

mu_o=4*pi*10^-7;
c=299792458;
eps_o=(1/(c*c*mu_o));
% i=[NPML_x:1:Nx-NPML_x];
% j=[NPML_y:1:Ny-NPML_y];


% f= (P*nab)E+(M*nab)H+(\ptP)X(mu_oH)-(\ptM)X(\eps_o E)
% fpx= (P*nab)E
%fJex=(\ptP)X(mu_oH)
%fJem=(\ptM)X(\eps_oE)
% fx=fpx+fJex+fJmx

% place fx ontop of Ex and Px i j will refer to i+1/2,j
Py_av=1/4*(Py(i,j)+Py(i+1,j)+Py(i,j-1)+Py(i+1,j-1));

fpx=Px(i,j).*((1/(2*dx)).*(Ex(i+1,j)-Ex(i-1,j)))+Py_av.*(1/(2*dy)).*(Ex(i,j+1)-Ex(i,j-1));

% fjx

Hz_av=1/4*(Hz(i,j)+Hz(i,j-1)+Hz_n_prev(i,j)+Hz_n_prev(i,j-1));
Jy_av=1/4*(Jey(i,j)+Jey(i+1,j)+Jey(i,j-1)+Jey(i+1,j-1));

fJex=Jy_av.*mu_o.*Hz_av;

% fjmx
Jmz_av=1/4*(Jmz(i,j)+Jmz(i,j-1)+Jmz_n_prev(i,j)+Jmz_n_prev(i,j-1));
Ey_avg=1/4*(Ey(i,j)+Ey(i+1,j)+Ey(i,j-1)+Ey(i+1,j-1));
fJmx=Jmz_av*eps_o.*Ey_avg;

fx(i,j)=fpx+fJex+fJmx;

end

