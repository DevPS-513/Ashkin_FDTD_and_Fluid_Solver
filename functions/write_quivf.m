function [ u_out,v_out ] = write_quivf( fname,x,y,u,v,m )

delete(fname)
fileID = fopen(fname,'w+');
[nrows]=length(x);




fprintf(fileID,' x \t y \t u \t v \t m \n ');

     for i=1:nrows
         
         check=1.*isnan([x(i) y(i) u(i) v(i) m(i)]);
         check=sum(check);
         if check==0
fprintf(fileID,' %.4f \t %.4f \t %.4f \t %.4f \t %.4f  \n ' ,x(i),y(i),u(i),v(i),m(i));
         end

     end


     

fclose(fileID)


end

