function [ Nf,xf,yf] = create_front( r,x,y,dx,dy )
% Define front based on density 
% This function will draw a single line countour around
% a given density profile

[X,Y]=meshgrid(x,y);

% Reduce step size for front
dx_f=dx/6;
dy_f=dy/6;

% Define start and end of [x,y]
x_start=x(1);
x_end=x(end);

y_start=y(1);
y_end=y(end);

xf_make=[x_start:dx_f:x_end];
yf_make=[y_start:dy_f:y_end];

[XF_MAKE, YF_MAKE]=meshgrid(xf_make,yf_make);

rho_f_make=interp2(X,Y,r,XF_MAKE,YF_MAKE);

front_make=contourc(rho_f_make,1);

xf_c=front_make(1,3:end-1)*dx_f;
yf_c=front_make(2,3:end-1)*dy_f;



Nf=length(xf_c); 
xf=zeros(1,Nf+2);
yf=zeros(1,Nf+2);

xf(1,2:end-1)=xf_c;
yf(1,2:end-1)=yf_c;



xf(1)=xf(Nf+1);
yf(1)=yf(Nf+1);
xf(Nf+2)=xf(2);
yf(Nf+2)=yf(2);


end

