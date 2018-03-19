
function [ er,g]=interface_create(x,y,cx1,avg,sig,h,er_slab,er_1)

% sig=1E-6;
% h=400e-9;
% avg=mean(y);

dx=x(2)-x(1);
dy=y(2)-y(1);
g=gauss_create(avg,sig,y);

g=g./max(g);

g=h.*g;
g=g+(cx1)-dx;

% figure
% plot(g,y,'color','black')
% hold on
% surf(x,y,er'-max(max(er)))
% shading flat
% view([ 0 90])re


[X,Y]=meshgrid(x,y);

Q=zeros(size(Y));

for i=1:length(x)
    
    for j=1:length(y)
     
        if((i-1)*dx>g(j))
            
          Q(j,i)=1;  
        end
      
        
        
        
        
    end
    
end

er=Q.*er_slab+1.*(Q~=1).*er_1;


% 
% figure(2)
% surf(Q)
% shading flat
% view([ 0 90])


end