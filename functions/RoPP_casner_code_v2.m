clear all
close all
clc
% include functions folder 
addpath functions

%% FIGURE INITIALIZATION

fig_count=10;
n_count=0;
s_d=get(0,'ScreenSize');  % Screen [0 0 width height]
sw=s_d(3);                % Screen width
sh=s_d(4);                % Screen height

fig_width = round(sw/3);
fig_height = round(sh/1.75);
long_fig_height = round(2*sh/3);
fig_x1 = 1;
fig_x2 = round(sw/2);
fig_y1 = round(sh/2.5);
fig_y2 = 1;

% Figure Positions
p1 = [fig_x1 fig_y1 fig_width fig_height];
p2 = [fig_x2 fig_y1 fig_width fig_height];


f_1=figure(1);
set(f_1,'name','Density')
set(f_1,'DefaulttextFontSize',14)
set(f_1,'OuterPosition',p1)

f_2=figure(2);
set(f_2,'name','Front')
set(f_2,'DefaulttextFontSize',14)
set(f_2,'OuterPosition',p2)



%% Domain Size 

Lx=1.0;
Ly=1.2;
Nx=42;
Ny=48;

% Grid step sizes

dx=Lx/Nx;
dy=Ly/Ny;
time_steps=10000; 

% Iteration parameters 
iterations=100;
tol=0.001;
beta=1.2;


%% Material and Geometry
rho1=1.0;
rho2=1.5;
min_rho=min([ rho1 rho2]);

mu1=.01;
mu2=2*mu1;
max_mu=max([mu1 mu2]);

% time step criterion

dt=(1/30)*(1/4)*min_rho*(min([ dx dy]))^2/(max_mu);
dt=.00125/2;
surf_tension=1;

% Initialize center of droplet
rad=.25;
xc=0.6;
yc=0.7; 

% Initialize velocity and pressure, and density
u=zeros(Nx+1,Ny+2);     %   u [m/s]
v=zeros(Nx+2,Ny+1);     %   v [m/s]
p=zeros(Nx+2,Ny+2);     %   p [pa]
u_star=zeros(Nx+1,Ny+2);    %   ut[m/s]
v_star=zeros(Nx+2,Ny+1);    %   vt[m/s]

grad_rx=zeros(Nx+2,Ny+2);   % density x-gradient
grad_ry=zeros(Nx+2,Ny+2);   % density y-gradient 

% Initialize Boundary conditions
u_wall_top=0;
u_wall_bottom=0;
v_wall_left=0;
v_wall_right=0;
% Initialize Gravity
gx=zeros(size(u));
gy=-100*ones(size(v));

% x and y axis, set to same as DNS solver to compare
x=zeros(1,Nx+2);
y=zeros(1,Ny+2);

for i=1:Nx+2;   
 x(i)=dx*(i-1.5);
end

for j=1:Ny+2;     
y(j)=dy*(j-1.5);
end


%% 2.0 Define Density

r=zeros(Nx+2,Ny+2)+rho1;
mu=zeros(Nx+2,Ny+2)+mu1;

% Rectangular indicies
x1=round(((xc-rad)-dx/2)/dx)+1;
x2=round(((xc+rad)-dx/2)/dx)+1;
y1=round(((yc-rad)-dy/2)/dy)+1;
y2=round(((yc+rad)-dy/2)/dy)+1;
 
% [Nr,r]=create_circle(r,xc,yc,rho2,r,dx,dy)

for i=2:Nx+1
    for j=2:Ny+1        

         if (( (x(i)-xc)^2+(y(j)-yc)^2)<rad^2)
             
            r(i,j)=rho2;
             mu(i,j)=mu2;
         end    
    end
end


%% CREATE r from FDTD
FDTD=load('N_data')
FDTD.er=FDTD.er';
[Nx_EM, Ny_EM]=size(FDTD.er);
dx_new=(x(end)-x(1))/(Nx_EM-1);
dy_new=(y(end)-y(1))/(Ny_EM-1);
i_sample=round((x+dx/2)/dx_new)+1;
j_sample=round((y+dy/2)/dy_new)+1;
er_mat=FDTD.er(i_sample,j_sample);

er_1=min(min(er_mat));
er_2=max(max(er_mat));

% er_mat(:,1)=er_1;
% er_mat(:,end)=er_1;
% er_mat(1,:)=er_1;
% er_mat(a

r=rho1.*(er_mat==er_1)+rho2.*(er_mat==er_2);
mu=mu1.*(er_mat==er_1)+mu2.*(er_mat==er_2);


r(:,1:2)=rho1;
r(:,end-1:end)=rho1;
r(1:2,:)=rho1;
r(end-1:end,:)=rho1;

%%  Initialize Front
   Nf=200;
%Create from density contour
 [Nf,xf,yf]=create_front_v2(r,x,y,dx,dy);
% %re-size front based on current grid
 [xf,yf,Nf]=resize_front(xf,yf,dx,dy,Nf);

%% END FDTD


        figure(1)
        surf(x,y,r');
        shading flat 
        xlabel('x-axis[m]'),ylabel('y-axis[m]')
        view([ 0 90])
        hold on
        plot3(xf,yf,(max(max(r))).*ones(size(xf)),'color','black')


% front velocities
    uf=zeros(1,Nf+2);
    vf=zeros(1,Nf+2);
% front tangent vectors     
    tx=zeros(1,Nf+2);
    ty=zeros(1,Nf+2);
% surface force
    fx=zeros(size(u));
    fy=zeros(size(v));   
    

amplitude=100;



%% 4.0 Time Loop
for n=1:time_steps
  
    %% Add external forc
    
    fx=0.*fx;
    fy=0.*fy;
    
    
%     [fx_ap,fy_ap]=ferro_force(x,y,dx,dy,amplitude);
%     
%     fx=fx+fx_ap(2:end,1:end);
%     fy=fy+fy_ap(1:end,2:end);
    
    
 % Refresh Boundary Conditions for wall velocities
 u(:,end)=2*(u_wall_bottom)-u(:,end-1);
 u(:,1)=2*(u_wall_top)-u(:,2); 
 v(1,:)=2*v_wall_left-v(2,:);
 v(end,:)=2*v_wall_right-v(end-1,:);   

% SURFACE TESNSION
    
   for q=1:Nf+1
       
ds=sqrt((xf(q+1)-xf(q))^2 +(yf(q+1)-yf(q))^2);
tx(q)=(xf(q+1)-xf(q))/ds;
ty(q)=(yf(q+1)-yf(q))/ds;

   end
tx(Nf+2)=tx(2);
ty(Nf+2)=ty(2);


% Disribute surface tension
 

    
for q=2:Nf+1
    
  
 % X-Surface Tension  
   % lands on u-velocites

iuf=floor(xf(q)/dx)+1;
juf=floor((yf(q)+dy/2)/dy)+1;

% Define interpolation square for u-locations
 
ax=xf(q)/dx+1-iuf;          % x distance to front point
ay=(yf(q)+dy/2)/dy+1-juf;   % y distance to front point


x1=0;
y1=0;
x2=1;
y2=1;
    
% Calculate all points   
    
    ntx=((tx(q)-tx(q-1))); % x normal vector from tangent vectors    
    surf_x_factor=surf_tension*(ntx)/(dx*dy);
    
% Bilinear weights

w_11=(x2-ax)*(y2-ay);
w_12=(x2-ax)*(ay-y1);
w_22=(ax-x1)*(ay-y1);
w_21=(ax-x1)*(y2-ay);
    
fx(iuf,juf)     =   fx(iuf,juf)     +w_11*surf_x_factor;
fx(iuf,juf+1)   =   fx(iuf,juf+1)   +w_12*surf_x_factor;
fx(iuf+1,juf+1) =   fx(iuf+1,juf+1) +w_22*surf_x_factor;
fx(iuf+1,juf)   =   fx(iuf+1,juf)   +w_21*surf_x_factor;



% Y-Surface Tension

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
 
    nty=(ty(q)-ty(q-1));  
    surf_y_factor=surf_tension*(nty)/(dx*dy);
    
    % Bring back weights
    
w_11=(x2-ax)*(y2-ay);
w_12=(x2-ax)*(ay-y1);
w_22=(ax-x1)*(ay-y1);
w_21=(ax-x1)*(y2-ay);

fy(ivf,jvf)     =fy(ivf,jvf)    +w_11*surf_y_factor;
fy(ivf,jvf+1)   =fy(ivf,jvf+1)  +w_12*surf_y_factor;
fy(ivf+1,jvf+1) =fy(ivf+1,jvf+1)+w_22*surf_y_factor;
fy(ivf+1,jvf)   =fy(ivf+1,jvf)  +w_21*surf_y_factor;


end
    
%% Start Navier-Stokes Algorithim
[Ax,Ay]=advect(u,v,dx,dy);

% Diffuse
[Dx,Dy]=diffuse_visc(u,v,dx,dy,mu);

% Intermediate velocities
[u_star,v_star]=intermediate_velocity(u,v,r,gx,gy,fx,fy,Ax,Ay,Dx,Dy,dt);

% pressure calculate
[p]=pressure_iterate_mat(u_star,v_star,dx,dy,dt,r,p,iterations,tol,beta);

% Correct the velocities
[u,v]=velocity_correct(u,v,u_star,v_star,dx,dy,dt,r,p);


%% Advect the front
% interpolate new front velocity
[uf,vf]=bilinear_interp_front(u,v,xf,yf,dx,dy,Nf);

% Update front location

for q=2:Nf+1
xf(q)=xf(q)+dt*uf(q);
yf(q)=yf(q)+dt*vf(q);

end

% Check front point spaceing
[xf,yf,Nf]=resize_front(xf,yf,dx,dy,Nf);

% Smooth Density gradient 
grad_rx=0.*grad_rx;
grad_ry=0.*grad_ry;

[grad_rx,grad_ry]=gradient_smoothing( rho1,rho2,grad_rx,grad_ry,xf,yf,dx,dy,Nf );

% Update Density

for k=1:iterations
      old_r=r;
for i=2:Nx+1       
    for j=2:Ny+1
        
        r_av=(1/4)*(r(i-1,j)+r(i+1,j)+r(i,j-1)+r(i,j+1));
        
        r(i,j)=r_av+(1/4)*dx*(grad_rx(i-1,j)-grad_rx(i,j))...
                   +(1/4)*dy*(grad_ry(i,j-1)-grad_ry(i,j));     
    end    
end
                if (max(max(abs(old_r-r)))<tol)
                    break
                end
end



%% Plot Density and front


if(n_count==fig_count)

    figure(2)
    plot(xf(1:Nf),yf(1:Nf))
    xlabel('x-axis[m]')
    ylabel('y-axis')
    xlim([x(1) x(end)])
    ylim([y(1) y(end)])
   
    
end

if(n_count==fig_count)
    n_count=0;
end
n

n_count=n_count+1;
end
