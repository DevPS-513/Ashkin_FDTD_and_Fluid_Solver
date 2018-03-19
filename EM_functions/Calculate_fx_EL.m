function [ fx] = Calculate_fx_EL( i,j,Hz,Hz_n_prev,Jmz,Jmz_n_prev_fr,Ex,Ey,Px,Py,Mz,Jey,Jex,dx,dy )

mu_o=4*pi*10^-7;
c=299792458;
eps_o=(1/(c*c*mu_o));
eta_o=sqrt(mu_o/eps_o);
fx=zeros(size(Hz));
% f_EL=(P*del)E+(M*del)H+JeXmu_oH+JmXeps_oE

% JeXmu_oH-JmXeps_oE

Hz_avg=(1/4)*(Hz(i,j)+Hz(i-1,j)+Hz_n_prev(i,j)+Hz_n_prev(i-1,j));
Jmz_avg=(1/4)*(Jmz(i,j)+Jmz(i-1,j)+Jmz_n_prev_fr(i,j)+Jmz_n_prev_fr(i-1,j));
fJx=mu_o.*Jey(i,j).*Hz_avg+eps_o.*Jmz_avg.*Ey(i,j);


% Px term

Px_avg=1/4*(Px(i,j)+Px(i-1,j)+Px(i,j+1)+Px(i-1,j+1));
Ex2=1/2*(Ex(i,j)+Ex(i,j+1));
Ex1=1/2*(Ex(i-1,j)+Ex(i-1,j+1));
fPx=Px_avg.*(1/dx).*(Ex2-Ex1);

% Py term

Ex2=1/2*(Ex(i,j+1)+Ex(i-1,j+1));
Ex1=1/2*(Ex(i,j)+Ex(i-1,j));
fPy=Py(i,j).*(1/dy).*(Ex2-Ex1);



fx(i,j)=fJx+fPx+fPy;

end

