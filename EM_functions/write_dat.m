function [  ] = write_dat( fname,x,y )

delete(fname)
fileID = fopen(fname,'w+');
[nrows] = length(x);
 for i = 1:nrows
     

         fprintf(fileID,'%d \t %d \n' ,x(i),y(i));
         

 end
fclose(fileID)


end

