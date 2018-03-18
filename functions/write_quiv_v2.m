function [  ] = write_quiv_v2( fname,u,v,tol,norm,skipx,skipy,edge_mat,max_find)

delete(fname)
fileID = fopen(fname,'w+');

% normalize
[Nx ,Ny]=size(u);


x=[0:1:Nx-1]./(Nx-1);
y=[0:1:Ny-1]./(Ny-1);




[nrows]=length(x);
[ncols]=length(y);

x1=edge_mat(1);
x2=edge_mat(2);
y1=edge_mat(3);
y2=edge_mat(4);

x1=round(x1.*nrows);
x2=round(x2.*nrows);
y1=round(y1.*ncols);
y2=round(y2.*ncols);


 max_f=max(max(sqrt(u(x1:x2,y1:y2).^2+v(x1:x2,y1:y2).^2)));
 
u=u./max_f;
v=v./max_f;

u=u.*norm;
v=v.*norm;
max_found=0;

fprintf(fileID,' x \t y \t u \t v \n ');
 for j = y1:skipy:y2
     for i=x1:skipx:x2
 
         check=max((sqrt(u(x1:x2,j).^2+v(x1:x2,j).^2)));
         if ( check>tol.*norm)
             
                    
             
fprintf(fileID,' %.4f \t %.4f \t %.4f \t %.4f \n ' ,x(i),y(j),u(i,j),v(i,j));
        end
         end
     end

% add row where maximum values is

if max_find==1

counter=1;
N_count=1;
 for j = y1:round(.8*nrows)
      check=sqrt(u(y1:y2,j).^2+v(y1:y2,j).^2);
      check=sum(1.*(check>=norm));

      
    
     for i=x1:skipx:x2
check2=sqrt(u(i,j).^2+v(i,j).^2);


         if ((check>=norm)&&(counter<=N_count)&&(check2>tol*norm))
             
fprintf(fileID,' %.4f \t %.4f \t %.4f \t %.4f \n ' ,x(i),y(j),u(i,j),v(i,j));
         end
         
         

     end
     
     
     
 end
end
 fclose(fileID)
end



     
 




