function [ xf,yf,Nf ] = resize_front(xf,yf,dx,dy,Nf)

xf_old=xf;
yf_old=yf;

c=1;
for q=2:Nf+1
  
    % q=2
    %c=1
    % if points within
    
   
ds=sqrt(  ((xf_old(q)-xf(c))/dx)^2+((yf_old(q)-yf(c))/dy)^2);
% if distance between new point at q, and olf point at c
if (ds > 0.5)
            c=c+1;

            xf(c)=0.5*(xf_old(q)+xf(c-1));
            yf(c)=0.5*(yf_old(q)+yf(c-1));
            
            c=c+1;
            
            xf(c)=xf_old(q);
            yf(c)=yf_old(q);    
elseif (ds<0.25) 
    
else
                c=c+1;
                xf(c)=xf_old(q);
                yf(c)=yf_old(q);
                
 end
 end
    
    

Nf=c-1;

xf(1)=xf(Nf+1);
yf(1)=yf(Nf+1);
xf(Nf+2)=xf(2);
yf(Nf+2)=yf(2);



end



