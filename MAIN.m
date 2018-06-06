clear all
close all
clc


% This code solves navier stokes equation for two fluids seperated by a
%  defined "front" it uses the front tracking tutorial found at https://www3.nd.edu/~gtryggva/MultiphaseDNS/DNS-Solver.pdf
%  and within the "Reference Documents" folder of this code.


% An example of how to plot the ouput data will be uploaded soon, in the
%  meantime you can load the output after the simulation via commands such
%  as "AB_down=load('./output_data/AB_down.mat')"
%  for the AB form and a beam pointed downward.
%  retreive simulation data as 
%  [dx dy Nx Ny dt Nt fig_count rho1 rho2 er_low er_high]=AB_down.parameters
%  then plot lets say the fy distribution at n=100 via
%  surf(AB_down.AB_fy_down_100)

% As well, after "beam_start", an FDTD code is run that solves for the
%  electromagnetic fields where the "front" between the two fluids is treated as
%  the interface between two media with two different refractive indicies.
%  the FDTD code has absorbing boundary conditions and only considers one
%  interface of the droplet. The force density is recalculated after
%  "beam_recalculate" time has passed



% Important points

% First and foremost, this is a continuum fluid model under the influence of
% various models for EM force density, but has not at this time been
% verified to represent the actual fluid velocity values that one would
% expect, and typically the simulation is unstable for anything not on the
% scale of meters, a CFD student may be able to read tryggvasons
% tutorial and then figure out why this is, but i'm not sure

% Make sure "beam_recalculate_time" is not too long. It should be 1 fluid
%   dt, but the simulation will take a drastically long time and you will
%   notice no appreciable change in the force density until the EM force
%   actually moves the front to a different shape, which may take several
%   miliseconds
% If gravity is turned on, the fluid will take the majority of the
%   simulation time just settling into a stationary state. It also needs this
%   time so that a flat surface is created. 
% The default setting is a [Lx X Ly] = [1m X 1m] liquid enviroment,
% the electromagnetic force is then scaled up to this size and normalized
% to a its peak value, which can then be multiplied

% Important parameters

% Nx and Ny are the x and y fluid resolution

% ratio_of_EM_force_to_gravity: The EM force is normalized, and 
% this constant determines how strong it should be considered in the fluid
% simulation relative to gravity. This is to make the deformation occur on
% meaningful time scales that the simulation can accomplish, and to make
% sure the force is not too small or too large.

addpath output_data
addpath EM_functions
addpath functions
addpath material_data


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
    seconds=1;
    femptoseconds=1e-15;
    mu_o=4*pi*10^-7;
    c=299792458;
    eps_o=(1/(c*c*mu_o));
    eta_o=sqrt(mu_o/eps_o);
    Kg=1;    
                   % The simulation is complete
    
%% FIGURE INITIALIZATION

    
 % Turn Figure plotting on or off    
    plot_on=1;                      % This will display Density over time
    plot_on_EM=0;                   % This will display the EM field
    save_mode=1;                    % This will save all workspace variables after
                 

    fig_count=1;
    n_count=0;
    s_d=get(0,'ScreenSize');  % Screen [0 0 width height] used to make figure sizes relative to the 
                              % screen size
    sw=s_d(3);                % Screen width
    sh=s_d(4);                % Screen height

% Figure Positions
% The figures are placed in a 3x3 array viewed on the screen

    




%% CFD Parameters

    Lx=2;
    Ly=2;
    Nx=72;
    Ny=72;
    source_direction=s_case;

    sim_time=3.125*seconds;
    dt=.0125/10;                  % CFD time step
    Nt=round(sim_time/dt);
    
    t=[0:1:Nt]*dt;

   
    beam_recalculate_time=0.2*seconds;          % Wait this amount of time before recalculateing the force distribution
    beam_recalculate=round(beam_recalculate_time/dt);
  
    beam_start_time=t(Nt)-8*beam_recalculate_time;  % this means for a total of 4 times, every .2 seconds
                                                    % the EM force desnity
                                                    % will be updated for
                                                    % the new geometry.
    beam_start=round(beam_start_time/dt);           % number of time steps before beam engages    
    
    beam_refresh_counter=0;                         % Counter that rolls over every
    
    
    % Grid step sizes
    dx=Lx/Nx;
    dy=Ly/Ny;

% Iteration parameters 

    iterations=50;
    tol=1;
    beta=1.2;
    


%% CFD DATA, DENSITY AND SURFACE TENSION

    rho1=1;                     % Density of fluid that represents air
    rho2=1.2;                   % Density of fluid that represents water 
        
    min_rho=min([ rho1 rho2]);
    rho_diff=abs(rho1-rho2)/2;

    surf_tension=0.65;            % Surface Tension Coefficient
    mu1=.04;
    mu2=2*mu1;
    max_mu=max([mu1 mu2]);




%% Initialize output matricies


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
    gy=-1*100*ones(size(v));

% x and y axis
    x=zeros(1,Nx+2);
    y=zeros(1,Ny+2);

for (i=1:Nx+2)  
    x(i)=dx*(i-1.5);
end

for (j=1:Ny+2)     
    y(j)=dy*(j-1.5);
end

%% FDTD and index of refraction data data

    er_high=1.33^2;         % permitivty where HIGH DENSITY is i.e water
    er_low=1.00^2;          % permitivity where LOW DENSITY is i.e air

    a=4000;                 % amplitude of the electric field
    LAMBDA=500*nm;          % wavelength of the source
    dx_EM=32*nm;    
    dy_EM=dx;
    beam_width=2*LAMBDA;
    x_size=4*beam_width;
    y_size=(Ly/Lx)*x_size;
    sim_cycles=1.0*y_size/(LAMBDA/max([er_low er_high]));
    
                            % this should only be after the fluid conforms
                            % to the container (i.e settles down)

    n_low=sqrt(er_high);
    n_high=sqrt(er_low);
    lam_low=LAMBDA/n_low;
    lam_high=LAMBDA/n_high;
    intercept=Ly/(1+lam_high/lam_low);  % Define where boundary between two liquids
                                        % should be

    ratio_of_EM_force_to_gravity=45;
% Initialize center of water source droplet droplet
    x_offset=7*dx;
    nx1=round(x_offset/dx)+1;
    nx2=round((x(end)-x_offset)/dx)+1;
    ny1=round(x_offset/dy)+1;
    ny2=round(intercept/dy)+1;

    xc=mean(x);
    yc=mean(y);
    rad=(x(end)-x(1))/5;
    r=rho1.*ones(Nx+2,Ny+2);
    
    % Find equivalent circle
    % pir^2=L*W,r=sqrt(L*w/pi)
        

        x_ellipse=round((round(Nx/2)-1));
        y_ellipse=round((ny2-ny1)/1.7);
        
        mid_y=round(2+y_ellipse);
        mid_x=round((Nx+2)/2);

        n_ellipse=3;
        
    for i=1:Nx+2
        for j=1:Ny+2
            
            if ( ( abs(((i-mid_x)/x_ellipse)^n_ellipse)+abs(((j-mid_y)/y_ellipse)^n_ellipse) )<1)
            r(i,j)=rho2;
            end
        end
    end
    
   % r(nx1:nx2,ny1:ny2)=rho2;

    [Nf,xf,yf]=create_front_v2(r,x,y,dx,dy);
 

    mu=mu1.*(r==rho1)+mu2.*(r==rho2);

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
           
    p1 = [1 sh/2 sw/3 sh/2];
    p2 = [sw/3 sh/2 sw/3 sh/2];
    p3 = [1 35 sw/3 sh/2];
    p4 = [sw/3 35 sw/3 sh/2];
    p5 = [sw/3+p3(3) sh/2 sw/3 sh/2];
      
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
    
    % Figure of fluid density
    figure(3)
    h_31=surf(x,y,r');
    xlim([-.1 1.1*Lx])
    ylim([-.1 1.1*Ly])
    caxis([ .9 1.3])
    shading flat
    view([ 0 90])
    hold on
    h_32=plot3(xf(3:Nf-4),yf(3:Nf-4),10.*ones(size(xf(3:Nf-4))),'color','black','linewidth',1.9);
    Nf
    length(xf)
    xlabel('x (m)')
    ylabel('y (m)')
    axis equal
    title(strcat(model,'_{',source_direction,'}'))
    
    
    
    
    
end

% Prepare figure variables 

    fig_counter=0;                      % Counter to set plot refresh
    N_fig=round(Nt/fig_count);  % number of steps until a figure refresh

%     fx_fig=zeros(Nx+1,Ny+2,N_fig);
%     fy_fig=zeros(Nx+2,Ny+1,N_fig);
    
%% 4.0 Time Loop

    fx_EM=0;
    fy_EM=0;
    eta_EM=0;
    Hz_EM=0;
    disp_max=0;   % Amount of displaecment in fluid domain to wait until
    xf_old=xf;
    yf_old=yf;

    [Nx_fx, Ny_fx]=size(u');
    [Nx_fy, Ny_fy]=size(v');

for n=1:Nt
    
    

% Erase Previous force values  
    fx=0.*fx;    
    fy=0.*fy;
 

% Start radiation beam if "beam_start" number of timsteps have passed
if n>beam_start
 
    xf_norm=(xf(1:Nf))./(x(end)-x(1));
    yf_norm=(yf(1:Nf))./(y(end)-y(1));
    c_e=0;  
  
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
    
    
    if (beam_refresh_counter>=beam_recalculate) 

        
        
       [fx_EM, fy_EM,eta_EM,Hz_EM,Sy_avg,fx_avg,fy_avg,x_EM_fx,y_EM_fx,x_EM_fy,y_EM_fy]=FDTD(Nx_fx,Ny_fx,Nx_fy,Ny_fy,x,y,xf_EM,yf_EM,er_low,er_high,LAMBDA,dx_EM,beam_width,x_size,y_size,a,sim_cycles,model,plot_on_EM,source_direction);
    beam_refresh_counter=0;

    % fy_EM and fx_EM are normalized to the maximum force value that
    % exists.
    % Scale force (note that gravity is '100') 
    fx_EM=fx_EM.*ratio_of_EM_force_to_gravity;
    fy_EM=fy_EM.*ratio_of_EM_force_to_gravity;
    end

end

% Reset disp max after it reaches one dex
    if (disp_max>(sqrt(dx^2+dy^2)/2)) 
        disp_max=0;
    end
    fx=fx+fx_EM;
    fy=fy+fy_EM;

    str = sprintf(['t(n) is \t' num2str(round(t(n),3)) '(s)\tand\t' num2str(round(100*n/Nt,2)) '\tperc. done']);
disp([str])

    beam_refresh_counter=beam_refresh_counter+1;
   
 % Refresh Boundary Conditions for wall velocities
%      u(:,end)=2*(u_wall_bottom)-u(:,end-1);
%      u(:,1)=2*(u_wall_top)-u(:,2); 
%      v(1,:)=2*v_wall_left-v(2,:);
%      v(end,:)=2*v_wall_right-v(end-1,:);   
     
     % tangential velocity at boundaries
    u(1:Nx+1,1)=        -u(1:Nx+1,2);
    u(1:Nx+1,Ny+2)=     -u(1:Nx+1,Ny+1);
    v(1,1:Ny+1)=        -v(2,1:Ny+1);
    v(Nx+2,1:Ny+1)=     -v(Nx+1,1:Ny+1);

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



% Check for maximum displacement the front has made, wait until it reaches
% one dx or dy before calling FDTD
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



%% Plot Density and front


if(n_count==fig_count)
    if plot_on==1
 % Figure 1-the front  
%    if n>beam_start
%        figure(1)
%        plot(xf_norm,yf_norm)
%        hold on
%        plot(xf_EM,yf_EM,'color','r')
%        hold on
%    scatter(x1_em,y1_em,'g');
%    hold on
%    scatter(x2_em,y2_em,'g');
%    xlim([ -.2 1.2]);
%    ylim([ -.2 1.2]);
% %    hold on
% %    quiver(x,y(2:end),zeros(size(fy_EM))',fy_EM'.*(dx)./max(max(fy_EM)));
%    hold off
%    end
 % Figure 3- density profile
 set(h_31,'ZData',r')
 set(h_32,'XData',xf(2:Nf),'YData',yf(2:Nf),'ZData',10.*ones(size(xf(2:Nf))))
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

data.(name_param)=[ dx dy Nx Ny dt Nt fig_count rho1 rho2 er_low er_high];
   
  n_count=0;
end
    
% n*dt

n_count=n_count+1;

clear xf_EM yf_EM; 
    
end

        if save_mode==1
            save(strcat(pwd,'./output_data/',model,'_',source_direction),'data');
            save(strcat(pwd,'./output_data/',model,'_',source_direction,'_workspace'));
        end

    end
end
    



