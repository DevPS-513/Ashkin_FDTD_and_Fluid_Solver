function [Vq]= bilinear_interp_mod(x,y,V,x_q,y_q,dx,dy)


Nx=length(x_q);
Ny=length(y_q);

Vq=zeros(Nx,Ny);
w_11=zeros(Nx,Ny);
w_12=zeros(Nx,Ny);
w_22=zeros(Nx,Ny);
w_21=zeros(Nx,Ny);


x_start=x(1);
y_start=y(1);




i=1:Nx;
    
    j=1:Ny;

        
        % find bottom left value on the fine grid
     xq=x_q(i);
     yq=y_q(j);
        
ip=floor((xq(i)-x_start)/dx)+1;
jp=floor((yq(j)-y_start)/dy)+1;
        

x1=x(ip);
x2=x1+dx;
y1=y(jp);
y2=y1+dy;



% for i=1:Nx
%     
%     for j=1:Ny


% Bilinear weights

w_11(i,j)=(x2'-xq')*(y2-yq)./(dx*dy);
 w_12(i,j)=(x2'-xq')*(yq-y1)./(dx*dy);
 w_22(i,j)=(xq'-x1')*(yq-y1)./(dx*dy);
w_21(i,j)=(xq'-x1')*(y2-yq)./(dx*dy);


%     end
%     
% end
% Solve for uf
%     Vq(i,j)=   w_11.*(V(ip,jp))      +w_12.*(V(ip,jp))...
%             +w_22.*(V(ip,jp))  +w_21.*(V(ip,jp));


          Vq(i,j)=   w_11(i,j).*(V(ip,jp))      +w_12(i,j).*(V(ip,jp+1))...
            +w_22(i,j).*(V(ip+1,jp+1))  +w_21(i,j).*(V(ip+1,jp));
 


end