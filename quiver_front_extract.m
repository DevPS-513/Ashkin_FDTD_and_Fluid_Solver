
function [ fxq,fyq,xfq,yfq ] = quiver_front_extract( xf,yf,fxf,fyf,xf1,xf2,yf1,Nf,x_sq)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



    % counter for new matrix
    c=1;
    for jf=2:Nf-1
        
        check=(xf(jf)<=xf2)&&(xf(jf)>=xf1)&&(xf(jf))&&(yf(jf)>yf1);
        
        if check==1
            fxq(c)=fxf(jf);
            fyq(c)=fyf(jf);
            xfq(c)=xf(jf);
            yfq(c)=yf(jf);
            c=c+1;
        end
        
    end

    fxq=interp1(xfq,fxq,x_sq);
    fyq=interp1(xfq,fyq,x_sq);
    yfq=interp1(xfq,yfq,x_sq);
    xfq=x_sq;
    
    % Lastly remove nans
    
   fxq(isnan(fxq))=0;
   fyq(isnan(fyq))=0;
end

