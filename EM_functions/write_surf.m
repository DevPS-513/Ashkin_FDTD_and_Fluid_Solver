function [  ] = write_surf( fname,x,y,u )

delete(fname)
fileID = fopen(fname,'w+');
[nrows]=length(x);
[ncols]=length(y);
 for j = 1:ncols
     for i=1:ncols
         
         
fprintf(fileID,' %.4f \t %.4f \t %.4f  \n ' ,x(i),y(j),u(i,j));

     end
     fprintf(fileID,'\n' );

     
 end
fclose(fileID)


end

