function [ er,x_size,y_size ] = convert_rho_to_eps(r,er_1,er_2,dx,dy)




[Nx, Ny]=size(r);
dx_new=(x(end)-x(1))/(Nx_EM-1);
dy_new=(y(end)-y(1))/(Ny_EM-1);
i_sample=round((x+dx/2)/dx_new)+1;
j_sample=round((y+dy/2)/dy_new)+1;
er_mat=FDTD.er(i_sample,j_sample);

er_1=min(min(er_mat));
er_2=max(max(er_mat));

% er_mat(:,1)=er_1;
% er_mat(:,end)=er_1;
% er_mat(1,:)=er_1;
% er_mat(a

r=rho1.*(er_mat==er_1)+rho2.*(er_mat==er_2);
mu=mu1.*(er_mat==er_1)+mu2.*(er_mat==er_2);




end

