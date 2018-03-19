function [uf,vf]= bilinear_interp_front(u,v,xf_mat,yf_mat,dx,dy,Nf)
uf=zeros(1,Nf);
vf=zeros(1,Nf);
for q=2:Nf+1

    
    % Define   
    xf=xf_mat(q);
    yf=yf_mat(q);
    % x&y-index
    

    
%% uf solve

% Find bottom left u-value
iuf=floor(xf/dx)+1;
juf=floor((yf+.5*dy)/dy)+1;

% define interpolaion box
x1=0;
y1=0;
x2=1;
y2=1;

ax=xf/dx+1-iuf;
ay=(yf+.5*dy)/dy+1-juf;

% Bilinear weights
w_11=(x2-ax)*(y2-ay);
w_12=(x2-ax)*(ay-y1);
w_22=(ax-x1)*(ay-y1);
w_21=(ax-x1)*(y2-ay);
% Solve for uf
    uf(q)=   w_11*(u(iuf,juf))      +w_12*(u(iuf,juf+1))...
            +w_22*(u(iuf+1,juf+1))  +w_21*(u(iuf+1,juf));

%% vf solve

% Index
ivf=floor((xf+.5*dx)/dx)+1;
jvf=floor(yf/dy)+1;

ax=(xf+.5*dx)/dx+1-ivf;
ay=yf/dy+1-jvf;


w_11=(x2-ax)*(y2-ay);
w_12=(x2-ax)*(ay-y1);
w_22=(ax-x1)*(ay-y1);
w_21=(ax-x1)*(y2-ay);
% Solve for uf
    vf(q)=   w_11*(v(ivf,jvf))      +w_12*(v(ivf,jvf+1))...
            +w_22*(v(ivf+1,jvf+1))  +w_21*(v(ivf+1,jvf));

        
        
end


end