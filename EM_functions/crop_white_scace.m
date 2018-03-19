function [ X ] = crop_white_space( X )
% input RGB image data,
% will output image data and crop outside white space

R=X(:,:,1);
G=X(:,:,2);
B=X(:,:,3);

W=1.*(R==0)+1.*(G==0)+1.*(B==0);

W=1.*(W~=0);

[row,col]=find(W>0);

crop_x1=min(col);
crop_x2=max(col);

crop_y1=min(row);
crop_y2=max(row);

X=X(crop_y1:crop_y2,crop_x1:crop_x2,:);
end

