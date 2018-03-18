function [ u_out,v_out ] = write_quiv( fname,x,y,u,v,m )

delete(fname)
fileID = fopen(fname,'w+');
[nrows]=length(x);
[ncols]=length(y);



fprintf(fileID,' x \t y \t u \t v \t m \n ');
 for j = 1:ncols
     for i=1:ncols
         
fprintf(fileID,' %.4f \t %.4f \t %.4f \t %.4f \t %.4f  \n ' ,x(i),y(j),u(i,j),v(i,j));

     end


     
 end
fclose(fileID)


end

