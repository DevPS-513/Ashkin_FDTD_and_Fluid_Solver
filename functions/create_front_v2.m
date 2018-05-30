function [ Nf,xf,yf] = create_front_v2( r,x,y,dx,dy )
% Define front based on density 
% This function will draw a single line countour around
% a given density profile


front_make=contourc(x,y,r',1);

xf_c=front_make(1,3:end-1);
yf_c=front_make(2,3:end-1);




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

