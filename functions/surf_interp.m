function [ h] = surf_interp( x,y,r,N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
r=r';
dx=x(2)-x(1);
dy=y(2)-y(1);
 
    x_new=[x(1)+dx/2:dx/N:x(end)-dx/2];
    y_new=[y(1)+dy/2:dy/N:y(end)-dy/2];
    
    [XN YN]=meshgrid(x_new,y_new);
    [XC YC]=meshgrid(x,y);    
    r_new=interp2(XC,YC,r,XN,YN,'spline');
    h=surf(x_new,y_new,r_new');
end

