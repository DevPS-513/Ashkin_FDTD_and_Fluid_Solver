
function [N,cx1,cx2,cy1,cy2]=create_multilayer(c_x,c_y,L_x,er_core,er_clad,N,dx,dy,N_core,d_core,d_clad)
%   input arguments
%   r-radius 
%   x- x cordinate of center, 
%   y-y coordinate of center
%   matrix to implement circle
%   eps_d permitivity of circle


%% Define L_y based on N_cell d_core and d_clad

L_y=(N_core-1)*d_core+(N_core)*d_clad-dy;

% create matrix at the center index of each core



r_nx=round(L_x/dx/2);
r_ny=round(L_y/dy/2);

m_y=round(c_y/dy);
m_x=round(c_x/dx);

% n as a function of r


cx1=m_x-r_nx;
cx2=m_x+r_nx;

cy1=m_y-r_ny;
cy2=m_y+r_ny;

n1=[cy1+round(d_clad/dy)+round(d_core/(2*dy))-1];
n2=[cy2-round(d_clad/dy)-round(d_core/(2*dy))+1];
    
n_jump= n1-cy1+2;
    
n_mat=[n1:n_jump:n2];

N(cx1:cx2,cy1:cy2)=er_clad;


for i=cx1:cx2    
    for j= (n1:n_jump:n2)    
        
    
        
        N(i,j-round(d_core/(2*dy))+1:j+round(d_core/(2*dy))-1)=er_core;
        

                  
             
        
    end   
end


end