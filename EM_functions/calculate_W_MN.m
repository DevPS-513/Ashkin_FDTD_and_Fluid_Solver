function [W]=calculate_W_MN(A,i,j,Dx,Dx_n_prev,Dy,Dy_n_prev,Bz,Bz_n_prev,dx,dy,dt,W)

% A is a a constant to multiply A.*(HXF)

% (Assume i,j refers to Bz locations)
% Average in time Dx and Dy values

Bz_x_2=1/4*(Bz(i+1,j)+Bz(i,j)+Bz_n_prev(i+1,j)+Bz_n_prev(i,j));

Bz_x_1=1/4*(Bz(i-1,j)+Bz(i,j)+Bz_n_prev(i-1,j)+Bz_n_prev(i,j));

Bz_y_2=1/4*(Bz(i,j+1)+Bz(i,j)+Bz_n_prev(i,j+1)+Bz_n_prev(i,j));

Bz_y_1=1/4*(Bz(i,j-1)+Bz(i,j)+Bz_n_prev(i,j-1)+Bz_n_prev(i,j));

X_term=(1/dx).*(  Dy(i+1,j).*Bz_x_2   -   Dy(i,j).*Bz_x_1 );
Y_term=(1/dy).*(  Dx(i,j+1).*Bz_y_2   -   Dx(i,j).*Bz_y_1 );

W(i,j)=W(i,j)-A.*dt.*(X_term-Y_term);


end

