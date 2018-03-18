function [ fx,fy ] = ferro_force( x,y,dx,dy,amp )
% clear all
% close all
% clc

% x=[-1:.03:-.1,.1:.03:1];
% y=x;


% extend the y-axis


% remember index for original axis

Nx_org=length(x);
Ny_org=length(y);



x_o=mean(x);
x=x-x_o;
Nx=length(x);
Ny=2*length(y);
y=[dy:Ny]*dy;
y_o=mean(y);
y=y-y_o;

ys=y(round(Ny/4));

y_o=mean(y);
x_o=mean(x);
sig_x=(1/25)*Nx*dx;
sig_y=(1/70)*Ny*dy;


 [X, Y]=meshgrid(y,x);




%% Create Gaussian


% gaussian parameters
x_term=(X-x_o).^2./(2*sig_x);
y_term=(Y-y_o).^2./(2*sig_y);

Z=amp.*exp(-1*(x_term+y_term));

theta=tan(X./Y);



Zx=Z.*cos(theta);
Zy=Z.*sin(theta);

%% Plot Gaussian
figure(2)
surf(x,y,Z')
shading flat
% caxis([ 0 100000])


% PLOT GRADIENT




i=1:Nx_org;
j=1:Ny_org;

fx=Zx(i,j);
fy=Zy(i,j);

figure(3)
quiver(fx,fy)

end

