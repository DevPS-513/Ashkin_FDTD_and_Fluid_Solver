
function [N2X,cx1,cx2,cy1,cy2]=create_box(L_x,L_y,c_x,c_y,eps_c,N2X,dx,dy)
%   input arguments
%   r-radius 
%   x- x cordinate of center, 
%   y-y coordinate of center
%   matrix to implement circle
%   eps_d permitivity of circle



r_nx=round(L_x/dx/2);
r_ny=round(L_y/dy/2);

m_y=round(c_y/dy);
m_x=round(c_x/dx);

% n as a function of r


cx1=m_x-r_nx;
cx2=m_x+r_nx;

cy1=m_y-r_ny;
cy2=m_y+r_ny;

for i=cx1:cx2    
    for j=cy1:cy2      
        
     
            
            
            N2X(i,j)=eps_c;            
             
        
    end   
end


end