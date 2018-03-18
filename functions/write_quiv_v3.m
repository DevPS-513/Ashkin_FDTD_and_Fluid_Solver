function [  ] = write_quiv_v3( fname,x,y,u,v,tol,norm,skipx,skipy,edge_mat,max_find,m)

delete(fname)
fileID = fopen(fname,'w+');

% normalize



% x=[0:1:Nx-1]./(Nx-1);
% y=[0:1:Ny-1]./(Ny-1);
% 
% [X,Y]=meshgrid(x,y);
% 
% 
% [nrows]=length(x);
% [ncols]=length(y);
% 
x1=1;
x2=length(x);
y1=1;
y2=length(y);
% 
% x1=round(x1.*nrows);
% x2=round(x2.*nrows);
% y1=round(y1.*ncols);
% y2=round(y2.*ncols);
% 
%  max_f=max(max(sqrt(u(x1:skipx:x2,y1:skipy:y2).^2+v(x1:skipx:x2,y1:skipy:y2).^2)));
%  
%  if max_f==0
%    u=u./1;
% v=v./1; 
%  else
%  
% u=u./max_f;
% v=v./max_f;
% 
%  end

u=u.*norm;
v=v.*norm;





fprintf(fileID,' x \t y \t u \t v \t m \n ');

if max_find==1
     for j = y1:1:y2
        for i=x1:1:x2
 
         check=max((sqrt(u(x1:x2,j).^2+v(x1:x2,j).^2)));
         
                 if ( (abs(check-norm)/norm)<=.001)
                fprintf(fileID,' %.4f \t %.4f \t %.4f \t %.4f \t %.4f\n ' ,x(i),y(j),u(i,j),v(i,j),m(i,j));
                 end
        
         
         end
     end    
end


 for j = y1:1:y2
     for i=x1:1:x2
 
         check=((sqrt(u(i,j).^2+v(i,j).^2)));
         
                if ( check>=tol)    
                fprintf(fileID,' %.4f \t %.4f \t %.4f \t %.4f \t %.4f\n ' ,x(i),y(j),u(i,j),v(i,j),m(i,j));
                end
     end
 end

% add row where maximum values is


 fclose(fileID)
end



     
 




