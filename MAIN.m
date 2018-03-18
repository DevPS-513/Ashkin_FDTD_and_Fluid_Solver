clear all
close all
clc

addpath output_data
addpath EM_functions
addpath functions
addpath material_data
addpath TE_EM_functions


models=cellstr(['AB ';'MN ';'EL ';'AMP';'Chu']);
s_cases=cellstr(['down';'up  ']);

mkdir('./output_data/')



for model_j=1:length(models)
    
    for case_j=1:length(s_cases)
   keep models s_cases model_j case_j
   clf
   close all
   
    model=char(models(model_j));
    s_case=char(s_cases(case_j));
%% SI units
    meters=1;
    nm=meters*1e-9;
    femptoseconds=1e-15;
    mu_o=4*pi*10^-7;
    c=299792458;
    eps_o=(1/(c*c*mu_o));
    eta_o=sqrt(mu_o/eps_o);
    Kg=1;    
    plot_on=1;
    plot_on_EM=0;
    save_mode=0;
%% FIGURE INITIALIZATION

    fig_count=2;
    n_count=0;
    s_d=get(0,'ScreenSize');  % Screen [0 0 width height]
    sw=s_d(3);                % Screen width
    sh=s_d(4);                % Screen height

    % Figure Positions
    p1 = [1 sh/2 sw/3 sh/2];
    p2 = [sw/3 sh/2 sw/3 sh/2];
    p3 = [1 35 sw/3 sh/2];
    p4 = [sw/3 35 sw/3 sh/2];
    p5 = [sw/3+p3(3) sh/2 sw/3 sh/2];
    p6 = [sw/3+p5(3) sh/2 sw/3 sh/2];
    if plot_on==1
    f_1=figure(1);
    set(f_1,'name','Front and <F>')
    set(f_1,'DefaulttextFontSize',14)
    set(f_1,'OuterPosition',p1)

    f_2=figure(2);
    set(f_2,'name','index in FDTD')
    set(f_2,'DefaulttextFontSize',14)
    set(f_2,'OuterPosition',p2)

    f_3=figure(3);
    set(f_3,'name','Density')
    set(f_3,'DefaulttextFontSize',14)
    set(f_3,'OuterPosition',p3)

    f_4=figure(4);
    set(f_4,'name','Hz_field')
    set(f_4,'DefaulttextFontSize',14)
    set(f_4,'OuterPosition',p4)

    f_5=figure(5);
    set(f_5,'name','\eta in FDTD')
    set(f_5,'DefaulttextFontSize',14)
    set(f_5,'OuterPosition',p5)

    f_6=figure(6);
    set(f_6,'name','Density')
    set(f_6,'DefaulttextFontSize',14)
    set(f_6,'OuterPosition',p6)

end



%% CFD Parameters
sim_case='EM';
mode_type='TE';
%sim_case='EM';
Lx=1;
Ly=1;
Nx=42;
Ny=42;
source_direction=s_case;

time_steps=500;
% Grid step sizes
dx=Lx/Nx;
dy=Ly/Ny;

% Iteration parameters 
    iterations=20;
    tol=.1;
    beta=1.2;


%% CFD DATA, DENSITY AND SURFACE TENSION
    rho1=1; rho2=1.2;
    min_rho=min([ rho1 rho2]);
    rho_diff=abs(rho1-rho2)/2;
    rho_water=1.2; % kg/m^3
    surf_tension=.7;
%surf_tension=1;
    mu1=.001;
    mu2=.001;
    max_mu=max([mu1 mu2]);
% time step criterion
    dt=(1/4)*min_rho*(min([ dx dy]))^2/(max_mu); % Minimum
    dt=.001;
    t=[0:1:time_steps]*dt;
    Nt=time_steps;


% Initialize velocity and pressure, and density
    F1x_avg=zeros(1,Nt);
    F1y_avg=zeros(1,Nt);
    F2x_avg=zeros(1,Nt);
    F2y_avg=zeros(1,Nt);
    u=zeros(Nx+1,Ny+2);         %   u [m/s]
    v=zeros(Nx+2,Ny+1);         %   v [m/s]
    p=zeros(Nx+2,Ny+2);         %   p [pa]
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
    grav=0;
    gy=-1*9.8*ones(size(v));

% x and y axis
    x=zeros(1,Nx+2);
    y=zeros(1,Ny+2);

    for i=1:Nx+2;   
     x(i)=dx*(i-1.5);
    end

    for j=1:Ny+2;     
    y(j)=dy*(j-1.5);
    end

%% FDTD data

    er_high=1.33^2;   % permitivty where HIGH DENSITY is
    er_low=1.0^2;       % permitivity where LOW DENSITY is

    a=-1;          % amplitude of the electric field
    f_scale=20;
    LAMBDA=500*nm;  % wavelength of the source
    dx_EM=18*nm;    
    dy_EM=dx;
    beam_width=2*LAMBDA;
    sos=0*LAMBDA;
    x_size=2.0*beam_width;
    y_size=(Ly/Lx)*x_size;
    sim_cycles=2.0*y_size/(LAMBDA/max([er_low er_high]));
    beam_start=1;

%% 2.0 Define Density and permitivity

    n_low=sqrt(er_high);
    n_high=sqrt(er_low);
    lam_low=LAMBDA/n_low;
    lam_high=LAMBDA/n_high;
    intercept=Ly/(1+lam_high/lam_low);
% intercept=intercept+.1;

% Initialize center of droplet
    x_offset=2*dx;
    nx1=round(x_offset/dx)+1;
    nx2=round((x(end)-x_offset)/dx)+1;
    ny1=round(2*dy/dy)+1;
    ny2=round(intercept/dy)+1;

    xc=mean(x);
    yc=mean(y);
    rad=(x(end)-x(1))/5;
    r=rho1.*ones(Nx+2,Ny+2);
    r(nx1:nx2,ny1:ny2)=rho2;

if strcmp(sim_case,'droplet')
    r=rho1.*ones(Nx+2,Ny+2);
    for i=1:Nx
        for j=1:Ny

            c1=(x(i)-xc)^2+(y(j)-yc)^2;

            if (c1<=rad^2)
               r(i,j)=rho2; 
            end

        end
    end

 y_bar_droplet=zeros(size(y));
 y_bar_drop=zeros(1,Nt);
 y_bar_droplet_theo=zeros(1,Nt);
 M_theo=pi*rad^2*rho_water; % assumeing length is 1 in z direction
    
end

    [Nf,xf,yf]=create_front_v2(r,x,y,dx,dy);
    mu=mu1.*(r==rho1)+mu2.*(r==rho2);

%%  Initialize Front
%Create from density contour
% %re-size front based on current grid
    [xf,yf,Nf]=resize_front(xf,yf,dx,dy,Nf);
% front velocities
    uf=zeros(1,Nf+2);
    vf=zeros(1,Nf+2);
% front tangent vectors     
    tx=zeros(1,Nf+2);
    ty=zeros(1,Nf+2);
% surface force
    fx=zeros(size(u));
    fy=zeros(size(v));     
 
    %% FIGURE CREATE
  
if plot_on==1
    % f_2=figure(1);
    % h_11=plot(xf(1:Nf)./(x(end)-x(1)),yf(1:Nf)./(x(end)-x(1)));
    % xlabel('x-axis[m]')
    % ylabel('y-axis')
    % xlim([x(1)-dx x(end)+dx])
    % ylim([y(1)-dy y(end)+dy])
    % line([ x(1) x(1)],[y(1) y(end)],'color','black')
    % line([ x(end) x(end)],[y(1) y(end)],'color','black')
    % line([ x(1) x(end)],[y(1) y(1)],'color','black')
    % line([ x(1) x(end)],[y(end) y(end)],'color','black')
    % axis equal
    % title(strcat(model,'_{',source_direction,'}'))
    figure(3)
    h_31=surf(x,y,r');
    xlim([-.1 1.1]*Lx)
    ylim([-.1 1.1]*Ly)
    caxis([ .9 1.3]*rho1)
    shading flat
    view([ 0 90])
    hold on
    h_32=plot3(xf(3:end-1),yf(3:end-1),10.*ones(size(xf(3:end-1))),'color','black','linewidth',1.9);
    h_yb=scatter3(x(round(Nx/2)),mean(y),10,'linewidth',2);
    h_yb_theo=scatter3(x(round(Nx/2)),mean(y),10,'linewidth',2,'MarkerFaceColor','black');

    xlabel('x')
    ylabel('y')
    axis equal
    colorbar east
    title(strcat(model,'_{',source_direction,'}'))

end

% store density,impedance,Hz field, fx and fy, for latert
    fig_counter=0;
    N_fig=round(time_steps/fig_count);
    r_fig=zeros(Nx+2,Ny+2,N_fig);
    u=zeros(Nx+1,Ny+2);     %   u [m/s]
    v=zeros(Nx+2,Ny+1);     %   v [m/s]
    fx_fig=zeros(Nx+1,Ny+2,N_fig);
    fy_fig=zeros(Nx+2,Ny+1,N_fig);




%% 4.0 Time Loop



    fx_EM=0;
    fy_EM=0;
    eta_EM=0;
    Ez_EM=0;
    disp_max=0;

    xf_old=xf;
    yf_old=yf;


for n=1:time_steps
    
% store old front
    
  
    fx=0.*fx;    
    fy=0.*fy;
 
if (( n>=beam_start))  
    
 
       %% convert front to a percentage of the grid
    xf_norm=(xf(1:Nf))./(x(end)-x(1));
    yf_norm=(yf(1:Nf))./(y(end)-y(1));
    c_e=0;
    
    % find front location at x=.2, x=.8
    

 
    for qf=2:Nf-1
    
         x1_em=.1;
    x2_em=.9;
    c1=(xf_norm(qf)>x1_em);
    c2=(xf_norm(qf)<x2_em);
    c3=(yf_norm(qf)>(.3*(ny2-ny1))/Ny);
    
    
    
    if(c1&&c2&&c3)
       c_e=c_e+1;
        xf_EM(c_e)=xf_norm(qf);
        yf_EM(c_e)=yf_norm(qf);   
        
        c_x1=((abs(xf_norm(qf)-x1_em))<(2/Nx));
        c_x2=((abs(xf_norm(qf)-x2_em))<(2/Nx));
        
        
        if c_x1
          y1_em=yf_norm(qf);  
        end
        
         if c_x2
          y2_em=yf_norm(qf);  
        end
       
    end
    
    % log y values at x1 and x2
    
    end
   
    
    xl_tack=[0:1/Nx:x1_em];
    xr_tack=[x2_em:1/Ny:1];
    yl_tack=y1_em.*ones(size(xl_tack));
    yr_tack=y2_em.*ones(size(xr_tack));

  xf_EM=[xl_tack,xf_EM,xr_tack];
  
  
  yf_EM=[yl_tack,yf_EM,yr_tack];
  
  

    
if (n==beam_start)
    er=er_low.*(r<(rho2-rho_diff))+er_high.*(r>=(rho2-rho_diff));
if strcmp(mode_type,'TE')
   
    [fx_EM, fy_EM,eta_EM,Hz_EM,Sy_avg,fx_avg,fy_avg]=casner_EM_code_interp_TE(x,y,er,er_low,er_high,LAMBDA,dx_EM,beam_width,x_size,y_size,a,sim_cycles,model,plot_on_EM,source_direction,sos);


    fx_EM(isnan(fx_EM))=0;
    fy_EM(isnan(fy_EM))=0;

    fx_EM=f_scale*fx_EM./max(max(sqrt(fx_EM.^2+fy_EM.^2)));
    fy_EM=f_scale*fy_EM./max(max(sqrt(fx_EM.^2+fy_EM.^2)));


    fx=fx+fx_EM(1:end-1,1:end);
    fy=fy+fy_EM(1:end,2:end);
end
   
%    fx_EM=0.*fx_EM;
%    fy_EM=0.*fy_EM;
   
disp_max/dx;

end

if (n>beam_start)
        fx_EM=0.*fx_EM;
        fy_EM=0.*fy_EM;
   
end




if (disp_max>(sqrt(dx^2+dy^2)/2)) 
    disp_max=0;
end

end
     
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
    
    xf_old(q)=xf(q);
    yf_old(q)=yf(q);

    xf(q)=xf(q)+dt*uf(q);
    yf(q)=yf(q)+dt*vf(q);

    disp_x(q)=abs(xf(q)-xf_old(q));
    disp_y(q)=abs(yf(q)-yf_old(q));
end
    xf(1)=xf(Nf+1);
    yf(1)=yf(Nf+1);
    xf(Nf+2)=xf(2);
    yf(Nf+2)=yf(2);


disp_max=disp_max+max([disp_x disp_y]);



% Check front point spaceinga
[xf,yf,Nf]=resize_front(xf,yf,dx,dy,Nf);

% Smooth Density gradient 
grad_rx=0.*grad_rx;
grad_ry=0.*grad_ry;

[grad_rx,grad_ry]=gradient_smoothing( rho2,rho1,grad_rx,grad_ry,xf,yf,dx,dy,Nf );

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

if ((n_count==fig_count)||n==1)
    if plot_on==1

 set(h_31,'ZData',r')
 set(h_32,'XData',xf(1:Nf),'YData',yf(1:Nf),'ZData',10.*ones(size(xf(1:Nf))))

% 
if n==1
figure(6)
hold off
surf(fy_EM')
shading flat
view([ 0 90])
hold on
caxis([min(min(fy_EM)) max(max(fy_EM))])
plot3(xf(1:Nf),yf(1:Nf),max(max(abs(fy))).*ones(size(xf(1:Nf))))
drawnow
% system('PAUSE')
end
drawnow
    end
    
 

    
name_xf=strcat(model,'_xf_',source_direction,'_',num2str(n));
name_yf=strcat(model,'_yf_',source_direction,'_',num2str(n));
name_fx=strcat(model,'_fx_',source_direction,'_',num2str(n));
name_fy=strcat(model,'_fy_',source_direction,'_',num2str(n));
name_fx_EM=strcat(model,'_fx_EM_',source_direction,'_',num2str(n));
name_fy_EM=strcat(model,'_fy_EM_',source_direction,'_',num2str(n));
name_eta_EM=strcat(model,'_eta_EM_',source_direction,'_',num2str(n));
name_Hz_EM=strcat(model,'_Hz_EM_',source_direction,'_',num2str(n));
name_Sy_av=strcat(model,'_Sy_svg_',source_direction,'_',num2str(n));
name_Nf=strcat(model,'_Nf_',source_direction,'_',num2str(n));
name_p=strcat(model,'_p_',source_direction,'_',num2str(n));
name_u=strcat(model,'_u_',source_direction,'_',num2str(n));
name_v=strcat(model,'_v_',source_direction,'_',num2str(n));



name_r=strcat(model,'_r_',source_direction,'_',num2str(n));
name_param='parameters';
data.(name_xf)=xf;
data.(name_yf)=yf;
data.(name_fx)=fx;
data.(name_fy)=fy;
data.(name_r)=r;
data.(name_fx_EM)=fx_EM;
data.(name_fy_EM)=fy_EM;
data.(name_eta_EM)=eta_EM;
data.(name_Hz_EM)=Hz_EM;
data.(name_Sy_av)=Sy_avg;
data.(name_Nf)=Nf;
data.(name_p)=p;
data.(name_u)=u;
data.(name_v)=v;



data.(name_param)=[ dx dy Nx Ny dt Nt fig_count rho1 rho2 er_low er_high];
   
  n_count=0;
end
    
% n*dt

n_count=n_count+1;

clear xf_EM yf_EM; 
    
end

if save_mode==1
    save(strcat(pwd,'./output_data/',mode_type,model,'_',source_direction),'data');
    save(strcat(pwd,'./output_data/',mode_type,model,'_',source_direction,'_workspace'));
end

    end

    end



