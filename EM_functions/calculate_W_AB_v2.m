function [W]=calculate_W_AB_v2(i,j,Ex,Ex_n_prev,Ey,Ey_n_prev,Hz,Hz_n_prev,dx,dy,dt,W)



% (Assume i,j refers to Hz locations)
% Average in time Ex and Ey values

Hz_x_2=1/4*(Hz(i+1,j)+Hz(i,j)+Hz_n_prev(i+1,j)+Hz_n_prev(i,j));

Hz_x_1=1/4*(Hz(i-1,j)+Hz(i,j)+Hz_n_prev(i-1,j)+Hz_n_prev(i,j));

Hz_y_2=1/4*(Hz(i,j+1)+Hz(i,j)+Hz_n_prev(i,j+1)+Hz_n_prev(i,j));

Hz_y_1=1/4*(Hz(i,j-1)+Hz(i,j)+Hz_n_prev(i,j-1)+Hz_n_prev(i,j));

X_term=(1/dx).*(  Ey(i+1,j).*Hz_x_2   -   Ey(i,j).*Hz_x_1 );
Y_term=(1/dy).*(  Ex(i,j+1).*Hz_y_2   -   Ex(i,j).*Hz_y_1 );

W(i,j)=W(i,j)-dt.*(X_term-Y_term);


end

