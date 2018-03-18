clear all
close all
clc
addpath EM_functions



dy=0.05;
    dx=dy;
x=[ -15:dx:15 ];
y=[-30:dy:30];



[X,Y] = meshgrid(x,y);

Nx=length(x);
Ny=length(y);

V=zeros(Nx,Ny);

for i=1:Nx
    
    for j=1:Ny
        
      V(i,j)=x(i)^2+(y(j)-6)^2 ; 
    end
    
end

% 
figure(1)
surf(x,y,V')
title('Original Sampling');
shading flat
view([ -90 0])
xq=[ -15:.05:15];
    yq=[-30:.05:30];

[Xq,Yq] = meshgrid(xq,yq);

for n=1:1
Vq = bilinear_interp_mod(x,y,V,xq(2:end-1),yq(2:end-1),dx,dy);
% Vq = bilinear_interp_mod(x,y,V,xq(2:end-1),yq(2:end-1),dx,dy);
% Vq = bilinear_interp_mod(x,y,V,xq(2:end-1),yq(2:end-1),dx,dy);
%Vq=interp2(X,Y,V',Xq,Yq);
figure(2)
hold on
surf(xq(2:end-1),yq(2:end-1),Vq');
title('Linear Interpolation Using courser Grid')
shading flat
view([ -90 0])


n
end