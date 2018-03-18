function [ fx,fy ] = tension_smoothing( rho1,rho2,fx,fy,xf,yf,dx,dy,Nf )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

for q=2:Nf+1
    
  
 %% X-GRADIENT SMOOTHING   
   % lands on u-velocites
    iuf=floor(xf(q)/dx)+1;
    juf=floor((yf(q)+dy/2)/dy)+1;
    
    % Define x-y square with v for fx
 
    ax=xf(q)/dx+1-iuf;
    ay=(yf(q)+dy/2)/dy+1-juf;
    
 
x1=0;
y1=0;
x2=1;
y2=1;
    
    % Calculate all points
    xf_p=xf(q+1);
    xf_m=xf(q-1);
    
    yf_p=yf(q+1);
    yf_m=yf(q-1);
    
    ds=sqrt((yf_p-yf_m)^2+(xf_p-xf_m)^2);    
    delta_y=yf_p-yf_m;    
    n21_x=-(0.5)*(delta_y)/ds;   
    
    grad_x_factor=(rho2-rho1)*ds*n21_x;
    
    % Bring back weights
    
w_11=(x2-ax)*(y2-ay)/(dx*dy);
w_12=(x2-ax)*(ay-y1)/(dx*dy);
w_22=(ax-x1)*(ay-y1)/(dx*dy);
w_21=(ax-x1)*(y2-ay)/(dx*dy);
    
grad_rx(iuf,juf)=grad_rx(iuf,juf)+w_11*grad_x_factor;
grad_rx(iuf,juf+1)=grad_rx(iuf,juf+1)+w_12*grad_x_factor;
grad_rx(iuf+1,juf+1)=grad_rx(iuf+1,juf+1)+w_22*grad_x_factor;
grad_rx(iuf+1,juf)=grad_rx(iuf+1,juf)+w_21*grad_x_factor;



%% Y-GRADIENT SMOOTHING



      % Define bottom left v-index
    ivf=floor((xf(q)+dx/2)/dx)+1;
    jvf=floor(yf(q)/dy)+1;
    
    % Define x-y square with v for fx
 
    ax=(xf(q)+dx/2)/dx+1-ivf;
    ay=(yf(q))/dy+1-jvf;

x1=0;
y1=0;
x2=1;
y2=1;
    
    % Calculate all points
    xf_p=xf(q+1);
    xf_m=xf(q-1);
    
    yf_p=yf(q+1);
    yf_m=yf(q-1);
    
    ds=sqrt((yf_p-yf_m)^2+(xf_p-xf_m)^2);

    delta_x=xf_p-xf_m;
    
   
    n21_y=0.5*(delta_x)/ds;
    
     grad_y_factor=(rho2-rho1)*ds*n21_y;
    
    % Bring back weights
w_11=(x2-ax)*(y2-ay)/(dx*dy);
w_12=(x2-ax)*(ay-y1)/(dx*dy);
w_22=(ax-x1)*(ay-y1)/(dx*dy);
w_21=(ax-x1)*(y2-ay)/(dx*dy);
    
    
grad_ry(ivf,jvf)=grad_ry(ivf,jvf)+w_11*grad_y_factor;
grad_ry(ivf,jvf+1)=grad_ry(ivf,jvf+1)+w_12*grad_y_factor;
grad_ry(ivf+1,jvf+1)=grad_ry(ivf+1,jvf+1)+w_22*grad_y_factor;
grad_ry(ivf+1,jvf)=grad_ry(ivf+1,jvf)+w_21*grad_y_factor;

end

end

