
function [N2X,cx1,cx2,cy1,cy2]=create_circle(r,c_x,c_y,eps_c,N2X,dx2,dy2)
%   input arguments
%   r-radius 
%   x- x cordinate of center, 
%   y-y coordinate of center
%   matrix to implement circle
%   eps_d permitivity of circle


[Nx2,Ny2]=size(N2X);
r_n=round(r/dx2);

m_y=round(c_y/dy2);
m_x=round(c_x/dx2);

% n as a function of r


cx1=m_x-r_n;
cx2=m_x+r_n;

cy1=m_y-r_n;
cy2=m_y+r_n;

for i=cx1:cx2    
    for j=cy1:cy2      
        
        check= (m_x-i)^2+(m_y-j)^2;
        
        if check<r_n^2
            
            
            N2X(i,j)=eps_c;            
        end           
        
    end   
end


end