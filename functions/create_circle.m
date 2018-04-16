
function [N,f]=create_circle(r,c_x,c_y,eps_c,N,dx,dy)
%   input arguments
%   r-radius 
%   x- x cordinate of center, 
%   y-y coordinate of center
%   matrix to implement circle
%   eps_d permitivity of circle


[Nx,Ny]=size(N);

dr=2*dx+2*dy;
r_n=r/dx;
r_nr=round(r_n);



m_y=round(c_y/dy);
m_x=Nx-round(c_x/dx);

% n as a function of r


cx1=m_x-ceil(r_n);
cx2=m_x+ceil(r_n);

cy1=m_y-ceil(r_n);
cy2=m_y+ceil(r_n);
l=1;

for i=cx1:cx2    
    for j=cy1:cy2      
        
        check= (m_x-i)^2+(m_y-j)^2;
        
        if check<=r_n^2
            
            
            N(i,j)=eps_c;            
        end
        
        if ( (sqrt(check)<=r_nr+.5*sqrt(2))&&(sqrt(check)>=r_nr-.5*sqrt(2)))
            
            f(l,1)=i;
            f(l,2)=j;
            l=l+1;
        end
        
    end   
end


end