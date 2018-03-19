function [ image_mat ] = center_y_at(y, image_mat,c_y )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[Nx,Ny]=size(image_mat);


[ c_3 ]=contourc(image_mat',1);
% findverticies from countour
[x1, I1]=min(c_3(1,2:end));
y1=c_3(2,I1+1);

[y2 I2]=max(c_3(2,2:end));
x2=c_3(1,I2+1);

[x3 I3]=max(c_3(1,2:end));
y3=c_3(2,I3+1);

[y4 I4]=min(c_3(2,2:end));
x4=c_3(1,I4+1);
% y2 is maximum y index
% y4 is minimum y index

remain_length=Ny-((y2+1)-(y4-1)+1);
half_length=round(remain_length/2);

object=image_mat(:,y4:y2);
[Nxo,Nyo]=size(object);

h1=half_length;
h2=half_length;

if (2*half_length+Nyo<Ny)
    h1=h1+1;
    
end

if (2*half_length+Nyo>Ny)
    h1=h1-1;
    
end


if c_y==mean(y)
image_mat=zeros(Nx,Ny);
image_mat(:,h1:h1+Nyo-1)=object;
else
    dy=abs(y(4)-y(3));
    offset=round((c_y-mean(y))/dy);
    image_mat=zeros(Nx,Ny);
image_mat(:,h1+offset:h1+offset+Nyo-1)=object;
end

