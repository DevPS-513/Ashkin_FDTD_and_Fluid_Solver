function [  ] = write_mat( fname,x,y,z )

delete(fname)
fileID = fopen(fname,'w+');
[nrows]=length(x);
[ncols]=length(y);
 for i = 1:nrows
     for j=1:ncols
         
         
fprintf(fileID,' %6.2f \t %6.2f \t %6.2f \n ' ,x(i),y(j),z(i,j));

     end
     fprintf(fileID,'\n' );

     
 end
fclose(fileID)


end

