function [  ] = saveas_cropped( h,name,ext )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

     saveas(h,name,ext);
         X=imread(strcat(name,'.',ext));
X=crop_white_space_v2(X,1,0,0);
imwrite(X,strcat(name,'.',ext));
clear X

end

