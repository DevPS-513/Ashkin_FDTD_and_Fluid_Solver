


function [fx_EM,fy_EM,eta_EM,Ez_EM,Sy_avg,fx_avg,fy_avg]=casner_EM_code_interp_TE(x_CFD,y_CFD,r,er_low,er_high,LAMBDA,dx,beam_width,x_size,y_size,a,sim_cycles,model,plot_on,source_direction,sos)


meters=1;
nm=meters*1e-9;
femptoseconds=1e-15;
mu_o=4*pi*10^-7;
c=299792458;
eps_o=(1/(c*c*mu_o));
eta_o=sqrt(mu_o/eps_o);
Kg=1;

% FIGURE INITIALIZATION

fig_count=150;              
n_count=0;


%% Simulation Parameters

f=c./LAMBDA;                    % Center frequency of pulse
T=1./f;                         % Center period

dy=dx;
dt=dx/(2*c);
% TEMPORAL GRID

  
sim_time=(sim_cycles*(1/f));
Nt=round(sim_time/dt);

% GEOMETRY


  
    NPML=.7*LAMBDA;                    % PML size
   
    
    NPML_x=round(NPML/dx);              % size of PML along x    
    NPML_y=NPML_x;              % size of PML along y
    
   
    Nx=round(x_size/dx)+2*NPML_x;           % Grid x length
    Ny=round(y_size/dy)+2*NPML_y;           % Grid y length
    
  
    x=[0:1:Nx-1]*dx;                        % x axis 
    y=[0:1:Ny-1]*dy;                        % y axis 
    t=[0:1:Nt-1]*dt;                      % time axis in [s]    
    
if strcmp(source_direction,'down')
source1_y=Ny-NPML_y-2-round(sos/dy);            % y-location of source
a_o=1*1.33^2;
S_yloc=source1_y-2;

end

if strcmp(source_direction,'up')
    source1_y=NPML_y+2+round(sos/dy);
    a_o=-1;
    S_yloc=source1_y+2;    
end

%% convert course grid to fine grid


% er_CFD=er_CFD(3:end-3,3:end-3);




%% convert CFD front percentages to EM front percentages, normalized to
% PML
nx1=NPML_x;
nx2=Nx-NPML_x+1;
ny1=NPML_y;
ny2=Ny-NPML_y+1;

% create axis on EM grid 
x_norm=x(nx1:nx2)-x(nx1);
x_norm=x_norm./(x(nx2)-x(nx1));

y_norm=y(ny1:ny2)-y(ny1);
y_norm=y_norm./(y(ny2)-y(ny1));

[Nx_CFD Ny_CFD]=size(r);

x_cfd=[0:1:Nx_CFD-1]./(Nx_CFD-1);
y_cfd=[0:1:Ny_CFD-1]./(Ny_CFD-1);
dx_cfd=1./(Nx_CFD-1);
dy_cfd=1./(Ny_CFD-1);


r_pad= zeros(Nx_CFD+2,Ny_CFD+2);
r_pad(2:end-1,2:end-1)=r;
r_pad(:,1)=er_low;
r_pad(:,end)=er_low;
r_pad(1,:)=er_low;
r_pad(end,:)=er_low;

x_cfd_pad=[-dx_cfd/2 x_cfd x_cfd(end)+dx_cfd/2];
y_cfd_pad=[-dy_cfd/2 y_cfd y_cfd(end)+dy_cfd/2];


[X_CFD_pad,Y_CFD_pad]=meshgrid(  x_cfd_pad,y_cfd_pad);
[XN YN]=meshgrid(x_norm,y_norm);

% buff r

er_interp=interp2(X_CFD_pad,Y_CFD_pad,r_pad,XN,YN);

er=er_low.*ones(Nx,Ny);
er(nx1:nx2,ny1:ny2)=er_interp;
 




%er_interp(isnan(er_interp))=er_1;




mr=ones(size(er));



    % Initialize drude matricies
    gamma_e=zeros(size(er));
    omega_e=zeros(size(er));
    gamma_m=zeros(size(er));
    omega_m=zeros(size(er));

    sigma_m=zeros(size(er));
    sigma_e=zeros(size(er));
    


%% PML Implement

    sigma_PML_e=0;                          % Taper PML conductivity
    sigma_PML_max=90000;  
    n_os=round(Nx/5);

    eps_1_left=er(NPML_x+n_os,:);        	% Sample Surrounding Medium, call this medium 1     
    mu_1_left=mr(NPML_x+n_os,:);
    % fill gap


    eps_1_right=er(Nx-NPML_x+1-n_os,:);        	% Sample Surrounding Medium, call this medium 1     
    mu_1_right=mr(Nx-NPML_x+1-n_os,:);
    % fill gap

     % fill gap
    eps_1_bottom=er(:,NPML_y+n_os);        	% Sample Surrounding Medium, call this medium 1     
    mu_1_bottom=mr(:,NPML_y+n_os,:);
    % fill gap


    
    
        eps_1_top=er(:,Ny-NPML_y+1-n_os);        	% Sample Surrounding Medium, call this medium 1     
    mu_1_top=mr(:,Ny-NPML_y+1-n_os);
     
    
    eps_2_left=(eps_1_left);                          % Solve for Medium 2 in PML, assuming eps2=eps1
    mu_2_left=(eps_2_left./eps_1_left).*mu_1_left;                % This will only work if mr=1 for surrounding materials
    
    eps_2_right=(eps_1_right);                          % Solve for Medium 2 in PML, assuming eps2=eps1
    mu_2_right=(eps_2_right./eps_1_right).*mu_1_right;                % This will only work if mr=1 for surrounding materials

    
       eps_2_top=(eps_1_top);                          % Solve for Medium 2 in PML, assuming eps2=eps1
    mu_2_top=(eps_2_top./eps_1_top).*mu_1_top;                % This will only work if mr=1 for surrounding materials
    
    eps_2_bottom=(eps_1_bottom);                          % Solve for Medium 2 in PML, assuming eps2=eps1
    mu_2_bottom=(eps_2_bottom./eps_1_bottom).*mu_1_bottom;                % This will only work if mr=1 for surrounding materials

    
    
    m_fac_left=(mu_2_left./(eps_2_left)).*(mu_o/eps_o);      % sigma_m=m_fac*sigma_e
    m_fac_right=(mu_2_right./(eps_2_right)).*(mu_o/eps_o);      % sigma_m=m_fac*sigma_e
    
      m_fac_top=(mu_2_top./(eps_2_top)).*(mu_o/eps_o);      % sigma_m=m_fac*sigma_e
    m_fac_bottom=(mu_2_bottom./(eps_2_bottom)).*(mu_o/eps_o);      % sigma_m=m_fac*sigma_e

    

    for j_d=1:Ny
    er(1:NPML_x+1+n_os,j_d)=er(NPML_x+1+n_os,j_d);
    er((Nx-NPML_x-1-(n_os)):Nx,j_d)=er(Nx-NPML_x-1-n_os,j_d);
    mr(1:NPML_x+1+n_os,j_d)=mr(NPML_x+1+n_os,j_d);
    mr((Nx-NPML_x-1-(n_os)):Nx,j_d)=mr(Nx-NPML_x-1-n_os,j_d);
        
    end
    
  
    
    for i_d=1:Nx
    er(i_d,1:NPML_y+1+n_os)=er(i_d,NPML_y+1+n_os);
    er(i_d,(Ny-NPML_y-(1+n_os)):Ny)=er(i_d,Ny-NPML_y-1-n_os);
    mr(i_d,1:NPML_y+1+n_os)=mr(i_d,NPML_y+1+n_os);
    mr(i_d,(Ny-NPML_y-(1+n_os)):Ny)=mr(i_d,Ny-NPML_y-1-n_os);  
    end
    
    
    slope_x=(sigma_PML_e-sigma_PML_max)/(NPML_x-1);
    slope_y=(sigma_PML_e-sigma_PML_max)/(NPML_y-1);
    b_x=sigma_PML_e-slope_x*NPML_x;
    b_y=sigma_PML_e-slope_y*NPML_y;
    
    
    % set x PML
for i=1:NPML_x
    
        sigma_e(i,:)=slope_x*i+b_x;        
        sigma_e(Nx-i+1,:)=slope_x*i+b_x;
        sigma_m(i,:)=m_fac_left.*sigma_e(i,:);
        sigma_m(Nx-i+1,:)=m_fac_right.*sigma_e(Nx-i+1,:);
    
end

    % set y PML
for j=1:NPML_y
    

        sigma_e(:,j)=slope_y*j+b_y+sigma_e(:,j);        
        sigma_e(:,Ny-j+1)=slope_y*j+b_y+sigma_e(:,Ny-j+1);
        sigma_m(:,j)=m_fac_bottom.*sigma_e(:,j);
        sigma_m(:,Ny-j+1)=m_fac_top.*sigma_e(:,Ny-j+1);
    
end
    
 %% Source

index_n=sqrt(er(round(Ny/2),source1_y)*mr(round(Ny/2),source1_y));
a=a./sqrt(index_n);
f_start=.2E15;                              % lowest frequency component   
f_end=1.2E15;                               % highest frequency component
FWHM=1/2.9e-15;                             % Full Width Half Maximum
N=50;                                       % Number of carrier signals
df=(f_end-f_start)/N;                       % Frequency spaceing
tau=5*T;                                    % shift pulse in time domain
G=gauss_pulse(tau,f,FWHM,N,f_start,f_end);  % Create Spectrum
norm_G=max((G(:,2)));                       % Normalize sourcespc/nm

G(:,2)=G(:,2)./norm_G;                      % Normalize..                            

% FDTD scaleing
t_avg=T;            
sig_t=.5*(1/f);   
gauss_t=1*exp(-1.*((t-t_avg).^2)/(2*sig_t^2));  

% Y-Source GAUSSIAN PROFILE
gauss_x_width=beam_width;
gauss_x_avg=mean(x);
sig_x=beam_width/1.5;
if strcmp(source_direction,'down')
sig_x=sig_x/(2);
end
if strcmp(source_direction,'up')
sig_x=1.33^2*sig_x/2;
end

N_x_gauss=round(gauss_x_width/dx);
x_gauss=[1:Nx]*dx;

gauss_x=gauss_create(gauss_x_avg,sig_x,x_gauss);
gauss_x=gauss_x./max(gauss_x);

% gauss_x_start=gauss_x_avg-1.5*sig_x;
% gauss_x_end=gauss_x_avg+1.5*sig_x;

gauss_x_start=x(NPML_x+1);
gauss_x_end=x(Nx-NPML_x);

nxg1=round(gauss_x_start/dx);
nxg2=round(gauss_x_end/dx);

%% OUTPUT Fields
i=4:Nx-4;
j=4:Ny-4;
i_W=i;
j_W=[4:1:source1_y-4, source1_y+4:1:Ny-4]; 

fx=zeros(Nx,Ny);
fy=zeros(Nx,Ny);

fy_along_center=zeros(1,Ny);
Fy_center=zeros(1,Nt);
Fy_system_avg=zeros(1,Nt);

eta_EM=sqrt(er(nx1:nx2,ny1:ny2)./mr(nx1:nx2,ny1:ny2))';


N_avg=round(5*T/dt);
avg_count=1;
fx_avg=zeros(Nx,Ny);
fy_avg=zeros(Nx,Ny);

fx_track=zeros(Nx,Ny,N_avg);
fy_track=zeros(Nx,Ny,N_avg);

Ez=zeros(Nx,Ny);
Pz=zeros(Nx,Ny);
Dz=zeros(Nx,Ny);
Dz_at_By=zeros(Nx,Ny);
Ez_at_By=zeros(Nx,Ny);

Hy=zeros(Nx,Ny);
By=zeros(Nx,Ny);
My=zeros(Nx,Ny);

Hx=zeros(Nx,Ny);
Bx=zeros(Nx,Ny);
Mx=zeros(Nx,Ny);

% dispersive current densities
Jz=zeros(Nx,Ny);
Jmy=zeros(Nx,Ny);
Jmx=zeros(Nx,Ny);



Fy_system_avg=zeros(1,Nt); 



% Pulse_Output_Variables
    W=zeros(Nx,Ny);                     %[J/m^3]        Energy Density
    M_pulse=zeros(1,Nt);                %[Kg]           Total pulse mass

    % pulse x-variables
    G_x=zeros(Nx,Ny);                    %[Kg/(m^2*s)]   Momentum Density
    G_x_n_prev=zeros(Nx,Ny);
    G_y_n_prev=zeros(Nx,Ny);


    G_y=zeros(Nx,Ny);                   %[Kg/(m^2*s)]   Momentum Density
  

 % Total variables
Tx=zeros(Nx,Ny);                                % Tx component
Ty=zeros(Nx,Ny);                                % Ty component

g_mech_x=zeros(Nx,Ny);                         %Momentum density
g_mech_y=zeros(Nx,Ny);                         %Momentum density


Dz_at_Bx=zeros(Nx,Ny);

%% Indexes

i=4:Nx-4;
j=4:Ny-4;

%% Hx Coeff

MX1=(1/dt-0.5*gamma_m(i,j))./(1/dt+0.5*gamma_m(i,j));
MX2=mu_o.*omega_m(i,j).^2./(2.*(1/dt+0.5.*gamma_m(i,j)));

X1=(mu_o.*mr(i,j)/dt-0.5*sigma_m(i,j)-MX2./2)./(mu_o.*mr(i,j)/dt+0.5*sigma_m(i,j)+MX2./2);
X2=1./(dy*((mu_o.*mr(i,j)/dt+0.5*sigma_m(i,j) +MX2./2)));
X3=(MX1+1)./(2.*(mu_o.*mr(i,j)/dt+0.5*sigma_m(i,j)+MX2./2));

BX1=mu_o.*mr(i,j);

%% Hy Coeff

MY1=(1/dt-0.5*gamma_m(i,j))./(1/dt+0.5*gamma_m(i,j));
MY2=mu_o*omega_m(i,j).^2./(2*(1/dt+0.5*gamma_m(i,j)));

Y1=(mu_o.*mr(i,j)/dt-0.5*sigma_m(i,j)-MY2/2)./(mu_o.*mr(i,j)/dt+0.5*sigma_m(i,j)+MY2./2);
Y2=1./(dx*((mu_o.*mr(i,j)/dt+0.5*sigma_m(i,j) +MY2./2)));
Y3=(MY1+1)./(2.*(mu_o.*mr(i,j)/dt+0.5*sigma_m(i,j)+MY2./2));

BY1=mu_o.*mr(i,j);

%% Ez Coeff


JZ1=( 1/dt-gamma_e(i,j)/2)./(1/dt+gamma_e(i,j)/2);
JZ2=eps_o*omega_e(i,j).^2./(2*(1/dt+gamma_e(i,j)/2));

Z1=(eps_o.*er(i,j)/dt-0.5*sigma_e(i,j)-JZ2./2)./(eps_o.*er(i,j)/dt+0.5*sigma_e(i,j)+JZ2./2);
Z2=1./(dx*((eps_o.*er(i,j)./dt+0.5*sigma_e(i,j)+JZ2./2)));
Z3=1./(dy*((eps_o.*er(i,j)./dt+0.5*sigma_e(i,j)+JZ2./2)));
Z4=(JZ1+1)./(2*((eps_o.*er(i,j)/dt+0.5*sigma_e(i,j)+JZ2./2)));

D1=eps_o.*er(i,j);
%% Indexes


% i=NPML_x+4:Nx-NPML_x-4;
% j=NPML_y+4:Ny-NPML_y-4;
% figure(4)
% xlabel(' Time ')
% ylabel(' Force(N)')
% hold on
% % legend('F_{tot}(t),F_{avg}(t)')
% 
% figure(4)
% fig_4_1=plot(t,F_system_y);
% hold on
% fig_4_2=plot(t,Fy_system_avg,'color','black');
% hold off
% xlabel('Time')
% ylabel(' F at center (N) ')
% 
% figure(3)
% subplot(1,2,1)
% fig_31=surf(x*1e6,y*1e6,Hz');
% % set(fig_3,'Colormap',[0.0431372560560703 0.517647087574005 0.780392169952393;0.0400560237467289 0.552100896835327 0.796078443527222;0.0369747914373875 0.586554646492004 0.811764717102051;0.033893559128046 0.621008396148682 0.82745099067688;0.0308123249560595 0.655462205410004 0.843137264251709;0.027731092646718 0.689916014671326 0.858823537826538;0.0246498603373766 0.724369764328003 0.874509811401367;0.0215686280280352 0.75882351398468 0.890196084976196;0.0184873957186937 0.793277323246002 0.905882358551025;0.0154061624780297 0.827731132507324 0.921568632125854;0.0123249301686883 0.862184882164001 0.937254905700684;0.00924369785934687 0.896638631820679 0.952941179275513;0.00616246508434415 0.931092441082001 0.968627452850342;0.00308123254217207 0.965546250343323 0.984313726425171;0 1 1;0.0555555559694767 1 1;0.111111111938953 1 1;0.16666667163372 1 1;0.222222223877907 1 1;0.277777791023254 1 1;0.333333343267441 1 1;0.388888895511627 1 1;0.444444447755814 1 1;0.5 1 1;0.555555582046509 1 1;0.611111104488373 1 1;0.666666686534882 1 1;0.722222208976746 1 1;0.777777791023254 1 1;0.833333313465118 1 1;0.888888895511627 1 1;0.944444417953491 1 1;1 1 1;1 0.923076927661896 0.923076927661896;1 0.846153855323792 0.846153855323792;1 0.769230782985687 0.769230782985687;1 0.692307710647583 0.692307710647583;1 0.615384638309479 0.615384638309479;1 0.538461565971375 0.538461565971375;1 0.461538463830948 0.461538463830948;1 0.384615391492844 0.384615391492844;1 0.307692319154739 0.307692319154739;1 0.230769231915474 0.230769231915474;1 0.15384615957737 0.15384615957737;1 0.0769230797886848 0.0769230797886848;1 0 0;0.972222208976746 0 0;0.944444417953491 0 0;0.916666686534882 0 0;0.888888895511627 0 0;0.861111104488373 0 0;0.833333313465118 0 0;0.805555582046509 0 0;0.777777791023254 0 0;0.75 0 0;0.722222208976746 0 0;0.694444417953491 0 0;0.666666686534882 0 0;0.638888895511627 0 0;0.611111104488373 0 0;0.583333313465118 0 0;0.555555582046509 0 0;0.527777791023254 0 0;0.5 0 0]);
% 
%        hold on
%        h_3=contour(x*1e6,y*1e6,er'-er_1,1,'color','red')
%        hold off
%        xlabel('x axis [{\mu}m]')
%        ylabel('y axis [{\mu}m]')
%      c_lim=    [-2*a/eta_o 2*a/eta_o];  
%        caxis(c_lim)
%        shading flat
% % Draw PML limes
%         line([NPML_x*dx  NPML_x*dx]*1e6, [(NPML_y)*dy (Ny-NPML_y)*dy]*1e6,'color','black');  
%         line([(Nx-NPML_x)*dx (Nx-NPML_x)*dx]*1e6, [(NPML_y)*dy (Ny-NPML_y)*dy]*1e6,'color','black');  
%         line([NPML_x*dx (Nx-NPML_x)*dx]*1e6, [(NPML_y)*dy (NPML_y)*dy]*1e6,'color','black');  
%         line([NPML_x*dx (Nx-NPML_x)*dx]*1e6, [(Ny-NPML_y)*dy (Ny-NPML_y)*dy]*1e6,'color','black'); %         
% % % Draw slab and pulse x- center of mass
% % if((x_bar_pulse(1,n)<(Nx-NPML_x)*dx)&&(y_bar_pulse(1,n)<(Ny-NPML_y)))

% view([ 0 90])
% title(' FDTD domain')
if plot_on==1
figure(2)
surf(([nx1:nx2]-nx1)./(nx2-nx1),([ny1:ny2]-ny1)./(ny2-ny1),er(nx1:nx2,ny1:ny2)')
xlim([-.2 1.2 ])
ylim([-.2 1.2])
view([ 0 90])
shading flat
title(' index n in FDTD')

hold off
figure(4)
h_41=surf(x,y,Ez');
shading flat
view([0 90])
hold off
figure(5)
surf(x*1e6,y*1e6,sqrt(mr./er)')
shading flat
view([0 90])
title(' \eta in FDTD')
line([NPML_x*dx  NPML_x*dx]*1e6, [(NPML_y)*dy (Ny-NPML_y)*dy]*1e6,[10 10],'color','black');  
line([(Nx-NPML_x)*dx (Nx-NPML_x)*dx]*1e6, [(NPML_y)*dy (Ny-NPML_y)*dy]*1e6,[10 10],'color','black');  
line([NPML_x*dx (Nx-NPML_x)*dx]*1e6, [(NPML_y)*dy (NPML_y)*dy]*1e6,[10 10],'color','black');  
line([NPML_x*dx (Nx-NPML_x)*dx]*1e6, [(Ny-NPML_y)*dy (Ny-NPML_y)*dy]*1e6,[10 10],'color','black'); % 
hold off


end

% 


eta_EM=sqrt(er(nx1:nx2,ny1:ny2)./mr(nx1:nx2,ny1:ny2))';

    for n= 1:Nt
      
          eta_r=sqrt(mr(round(Nx/2),source1_y)./er(round(Nx/2),source1_y));
        Ez(nxg1:nxg2,source1_y)=a_o*gauss_x(nxg1:nxg2)'.*(a/(2*eta_o)).*sin(2*pi*f*(t(n)))+ Ez(nxg1:nxg2,source1_y);            % X- Guassian %         
    
  
        
     Hx_n_prev=Hx;
        Hy_n_prev=Hy;
        Ez_n_prev=Ez;        
        
        Bx_n_prev=Bx;
        By_n_prev=By;
        Dz_n_prev=Dz;
        
        Jz_n_prev=Jz;
        Jmx_n_prev=Jmx;
        Jmy_n_prev=Jmy;    

%% H Update (n+1/2)
Hx(i,j)=X1.*Hx(i,j)-X2.*(Ez(i,j)-Ez(i,j-1))-X3.*Jmx(i,j);
Jmx(i,j)=MX1.*Jmx(i,j)+MX2.*(Hx(i,j)+Hx_n_prev(i,j));
Mx(i,j)=Mx(i,j)+(dt/2).*(Jmx(i,j)+Jmx_n_prev(i,j));
Bx(i,j)=BX1.*Hx(i,j)+Mx(i,j);

Hy(i,j)=Y1.*Hy(i,j)+Y2.*(Ez(i,j)-Ez(i-1,j))-Y3.*Jmy(i,j);        
Jmy(i,j)=MY1.*Jmy(i,j)+MY2.*(Hy(i,j)+Hy_n_prev(i,j));        
My(i,j)=My(i,j)+(dt/2).*(Jmy(i,j)+Jmy_n_prev(i,j));        
By(i,j)=BY1.*Hy(i,j)+My(i,j);     %Hz(source1_x,nyg1:nyg2)=gauss_y(nyg1:nyg2).*(a/(2*eta_o)).*cos(2*pi*f*(t(n)))+ Hz(source1_x,nyg1:nyg2);            % X- Guassian %        
 



    
%% MN (n+1/2) CALCULATIONS
if strcmp(model,'MN')
        [W]=calculate_W_MN(c^2,i,j,Bx,Bx_n_prev,By,By_n_prev,-1*Dz,-1*Dz_n_prev,dx,dy,dt,W);



% x-calculations
Dz_at_By(i,j)=.25*(Dz(i,j)+Dz(i-1,j)+Dz_n_prev(i,j)+Dz_n_prev(i-1,j));
[Tx,t1x,t2x,t3x,t4x,t5x]=calculate_Tx_MN_TE(i,j, Bx,By,Hx,Bx_n_prev,By_n_prev,Hx_n_prev,Hy_n_prev,Hy,Ez,Dz,dx,dy );
G_x_n_prev(i,j)=G_x(i,j);
G_x(i,j)=-1.*1.*Dz_at_By(i,j).*By(i,j);
% y-calculations
% y-calcualtions
Dz_at_Bx(i,j)=.25*(Dz(i,j)+Dz(i,j-1)+Dz_n_prev(i,j)+Dz_n_prev(i,j-1));
[Ty,t1y,t2y,t3y,t4y,t5y]=calculate_Ty_MN_TE(i,j, Bx,By,Hx,Bx_n_prev,By_n_prev,Hx_n_prev,Hy_n_prev,Hy,Ez,Dz,dx,dy );
G_y_n_prev(i,j)=G_y(i,j);
G_y(i,j)=Dz_at_Bx(i,j).*Bx(i,j);

%W(i,j)=W(i,j)-(c.^2.*(dt/dx)*(G_x(i+1,j)-G_x(i,j))-c.^2*(dt/dy).*(G_y(i,j+1)-G_y(i,j)));
      % POWER CALCULATIONS
    Sx(i,j)=(G_x(i,j)).*c^2;
    Sy(i,j)=(G_y(i,j)).*c^2;
    
% if n>round(4*T/dt)
% Sy_flux(1,n)=sum(Sy(:,S_yloc)*dy);
% Sy_avg(1,n)=sum(Sy_flux(n-round(4*T/dt):n))/(round(4*T/dt));
% end

if n>1
Sy_flux(1,n)=sum(Sy(:,S_yloc)*dy);
Sy_avg(1,n)=mean(Sy_flux(1:n));
end

end

if strcmp(model,'AB')
        [W]=calculate_W_AB_v2(i,j,Hx,Hx_n_prev,Hy,Hy_n_prev,-1*Ez,-1*Ez_n_prev,dx,dy,dt,W);      




% x-calculations
Ez_at_By(i,j)=.25*(Ez(i,j)+Ez(i-1,j)+Ez_n_prev(i,j)+Ez_n_prev(i-1,j));
[Tx,t1x,t2x,t3x,t4x,t5x]=calculate_Tx_AB_TE(i,j, Bx,By,Hx,Bx_n_prev,By_n_prev,Hx_n_prev,Hy_n_prev,Hy,Ez,Dz,dx,dy );
G_x_n_prev(i,j)=G_x(i,j);
G_x(i,j)=-1.*1.*Ez_at_By(i,j).*Hy(i,j)./c^2;  
% y-calcualtions
Ez_at_Bx(i,j)=.25*(Ez(i,j)+Ez(i,j-1)+Ez_n_prev(i,j)+Ez_n_prev(i,j-1));
[Ty,t1y,t2y,t3y,t4y,t5y]=calculate_Ty_AB_TE(i,j, Bx,By,Hx,Bx_n_prev,By_n_prev,Hx_n_prev,Hy_n_prev,Hy,Ez,Dz,dx,dy );
G_y_n_prev(i,j)=G_y(i,j);
G_y(i,j)=Ez_at_Bx(i,j).*Hx(i,j)./c^2;

      % POWER CALCULATIONS
    Sx(i,j)=(G_x(i,j)).*c^2;
    Sy(i,j)=(G_y(i,j)).*c^2;
    
% if n>round(4*T/dt)
% Sy_flux(1,n)=sum(Sy(:,S_yloc)*dy);
% Sy_avg(1,n)=sum(Sy_flux(n-round(4*T/dt):n))/(round(4*T/dt));
% end

if n>1
Sy_flux(1,n)=sum(Sy(:,S_yloc)*dy);
Sy_avg(1,n)=mean(Sy_flux(1:n));
end


%W(i_W,j)=W(i_W,j)-(c.^2.*(dt/dx)*(G_x(i_W+1,j)-G_x(i_W,j))-c.^2*(dt/dy).*(G_y(i_W,j+1)-G_y(i_W,j)));


end


if strcmp(model,'EL')
        [W]=calculate_W_AB_v2(i,j,Hx,Hx_n_prev,Hy,Hy_n_prev,-1*Ez,-1*Ez_n_prev,dx,dy,dt,W);




% x-calculations
Ez_at_By(i,j)=.25*(Ez(i,j)+Ez(i-1,j)+Ez_n_prev(i,j)+Ez_n_prev(i-1,j));
[Tx,t1x,t2x,t3x,t4x,t5x]=calculate_Tx_EL_TE(i,j, Bx,By,Hx,Bx_n_prev,By_n_prev,Hx_n_prev,Hy_n_prev,Hy,Ez,Dz,dx,dy );
G_x_n_prev(i,j)=G_x(i,j);
G_x(i,j)=-1.*1.*Ez_at_By(i,j).*Hy(i,j)./c^2;  
% y-calcualtions
Ez_at_Bx(i,j)=.25*(Ez(i,j)+Ez(i,j-1)+Ez_n_prev(i,j)+Ez_n_prev(i,j-1));
[Ty,t1y,t2y,t3y,t4y,t5y]=calculate_Ty_EL_TE(i,j, Bx,By,Hx,Bx_n_prev,By_n_prev,Hx_n_prev,Hy_n_prev,Hy,Ez,Dz,dx,dy );
G_y_n_prev(i,j)=G_y(i,j);
G_y(i,j)=Ez_at_Bx(i,j).*Hx(i,j)./c^2;

      % POWER CALCULATIONS
    Sx(i,j)=(G_x(i,j)).*c^2;
    Sy(i,j)=(G_y(i,j)).*c^2;
    
% if n>round(4*T/dt)
% Sy_flux(1,n)=sum(Sy(:,S_yloc)*dy);
% Sy_avg(1,n)=sum(Sy_flux(n-round(4*T/dt):n))/(round(4*T/dt));
% end

if n>1
Sy_flux(1,n)=sum(Sy(:,S_yloc)*dy);
Sy_avg(1,n)=mean(Sy_flux(1:n));
end


%W(i_W,j)=W(i_W,j)-(c.^2.*(dt/dx)*(G_x(i_W+1,j)-G_x(i_W,j))-c.^2*(dt/dy).*(G_y(i_W,j+1)-G_y(i_W,j)));


end

if strcmp(model,'Chu')
        [W]=calculate_W_AB_v2(i,j,Hx,Hx_n_prev,Hy,Hy_n_prev,-1*Ez,-1*Ez_n_prev,dx,dy,dt,W);


% x-calculations
Ez_at_By(i,j)=.25*(Ez(i,j)+Ez(i-1,j)+Ez_n_prev(i,j)+Ez_n_prev(i-1,j));
[Tx,t1x,t2x,t3x,t4x,t5x]=calculate_Tx_Chu_TE(i,j, Bx,By,Hx,Bx_n_prev,By_n_prev,Hx_n_prev,Hy_n_prev,Hy,Ez,Dz,dx,dy );
G_x_n_prev(i,j)=G_x(i,j);
G_x(i,j)=-1.*1.*Ez_at_By(i,j).*Hy(i,j)./c^2;  
% y-calcualtions
Ez_at_Bx(i,j)=.25*(Ez(i,j)+Ez(i,j-1)+Ez_n_prev(i,j)+Ez_n_prev(i,j-1));
[Ty,t1y,t2y,t3y,t4y,t5y]=calculate_Ty_Chu_TE(i,j, Bx,By,Hx,Bx_n_prev,By_n_prev,Hx_n_prev,Hy_n_prev,Hy,Ez,Dz,dx,dy );
G_y_n_prev(i,j)=G_y(i,j);
G_y(i,j)=Ez_at_Bx(i,j).*Hx(i,j)./c^2;

      % POWER CALCULATIONS
    Sx(i,j)=(G_x(i,j)).*c^2;
    Sy(i,j)=(G_y(i,j)).*c^2;
    
% if n>round(4*T/dt)
% Sy_flux(1,n)=sum(Sy(:,S_yloc)*dy);
% Sy_avg(1,n)=sum(Sy_flux(n-round(4*T/dt):n))/(round(4*T/dt));
% end

if n>1
Sy_flux(1,n)=sum(Sy(:,S_yloc)*dy);
Sy_avg(1,n)=mean(Sy_flux(1:n));
end


%W(i_W,j)=W(i_W,j)-(c.^2.*(dt/dx)*(G_x(i_W+1,j)-G_x(i_W,j))-c.^2*(dt/dy).*(G_y(i_W,j+1)-G_y(i_W,j)));


end


if strcmp(model,'AMP')
            [W]=calculate_W_MN(eps_o*c^2,i,j,Bx,Bx_n_prev,By,By_n_prev,-1*Ez,-1*Ez_n_prev,dx,dy,dt,W);    
   

              


% x-calculations
Ez_at_By(i,j)=.25*(Ez(i,j)+Ez(i-1,j)+Ez_n_prev(i,j)+Ez_n_prev(i-1,j));
[Tx,t1x,t2x,t3x,t4x,t5x]=calculate_Tx_AMP_TE(i,j, Bx,By,Hx,Bx_n_prev,By_n_prev,Hx_n_prev,Hy_n_prev,Hy,Ez,Dz,dx,dy );
G_x_n_prev(i,j)=G_x(i,j);
G_x(i,j)=-1.*eps_o.*Ez_at_By(i,j).*By(i,j);  
% y-calcualtions
Ez_at_Bx(i,j)=.25*(Ez(i,j)+Ez(i,j-1)+Ez_n_prev(i,j)+Ez_n_prev(i,j-1));
[Ty,t1y,t2y,t3y,t4y,t5y]=calculate_Ty_AMP_TE(i,j, Bx,By,Hx,Bx_n_prev,By_n_prev,Hx_n_prev,Hy_n_prev,Hy,Ez,Dz,dx,dy );
G_y_n_prev(i,j)=G_y(i,j);
G_y(i,j)=eps_o.*Ez_at_Bx(i,j).*Bx(i,j);

      % POWER CALCULATIONS
    Sx(i,j)=(G_x(i,j)).*c^2;
    Sy(i,j)=(G_y(i,j)).*c^2;
    
% if n>round(4*T/dt)
% Sy_flux(1,n)=sum(Sy(:,S_yloc)*dy);
% Sy_avg(1,n)=sum(Sy_flux(n-round(4*T/dt):n))/(round(4*T/dt));
% end

if n>1
Sy_flux(1,n)=sum(Sy(:,S_yloc)*dy);
Sy_avg(1,n)=mean(Sy_flux(1:n));
end


%W(i_W,j)=W(i_W,j)-(c.^2.*(dt/dx)*(G_x(i_W+1,j)-G_x(i_W,j))-c.^2*(dt/dy).*(G_y(i_W,j+1)-G_y(i_W,j)));


end

 %% Ez uPDATE  (N+1)               
Ez(i,j)=Z1.*Ez(i,j)+Z2.*(Hy(i+1,j)-Hy(i,j))-Z3.*(Hx(i,j+1)-Hx(i,j))-Z4.*Jz(i,j);
Jz(i,j)=JZ1.*Jz(i,j)+JZ2.*(Ez(i,j)+Ez_n_prev(i,j));
Pz(i,j)=Pz(i,j)+(dt./2).*(Jz(i,j)+Jz_n_prev(i,j));
Dz(i,j)=D1.*Ez(i,j)+Pz(i,j);


       % Ex(nxg1:nxg2,source1_y)=a_o*gauss_x(nxg1:nxg2)'.*(a/(2)).*sin(2*pi*f*(t(n)))+ Ex(nxg1:nxg2,source1_y);            % X- Guassian %         


    y_bar_pulse(1,n)=(1/M_pulse(1,n)).*sum(y(j).*sum((W(i,j)./c^2).*dx)).*dy;
%% Slab Calculations

%     fx_T(i,j)=-Tx(i,j)-1*(1/dt)*(G_x(i,j)-G_x_n_prev(i,j));
%     fy_T(i,j)=-Ty(i,j)-1*(1/dt)*(G_y(i,j)-G_y_n_prev(i,j));
    g_mech_x_prev=g_mech_x;
g_mech_y_prev=g_mech_y;


    g_mech_x(i,j)=g_mech_x(i,j)-dt.*Tx(i,j)-1*(G_x(i,j)-G_x_n_prev(i,j));
    g_mech_y(i,j)=g_mech_y(i,j)-dt.*Ty(i,j)-1*(G_y(i,j)-G_y_n_prev(i,j));
    

    
    
 
 
 
fx(i,j)=(g_mech_x(i,j)-g_mech_x_prev(i,j))/dt;
fy(i,j)=(g_mech_y(i,j)-g_mech_y_prev(i,j))/dt;

%     fx(:,source1_y-6:source1_y+6)=0;
%     fy(:,source1_y-6:source1_y+6)=0;
    
   
%     
%     F_system_y(1,n)=sum(sum(g_mech_y(i,j)-g_mech_y_prev(i,j)))*dx*dy;
%     Fy_system_avg(1,n)=sum(F_system_y(1:n))/n;

     F_system_y(1,n)=sum(fy(round(Nx/2),j))*dy;
     
     
   
          if avg_count<=N_avg
         fx_track(:,:,avg_count)=fx;
         fy_track(:,:,avg_count)=fy;
         avg_count=avg_count+1;
         if avg_count==N_avg
             
            avg_count=1; 
         end
     end
if (n>(round(5*T/dt)+1))
    % Fy_system_avg(1,n)=mean(F_system_y(1,n-round(5*T/dt):n));
   
%Fy_system_avg(1,n)=mean(sum(sum((fy_track(round(Nx/2),j,:))*dx*dy)));    
Fy_system_avg(1,n)=mean(F_system_y(1,n-round(5*T/dt):n));

%fy_avg=fy_avg./max(max(fy_avg));

 
     
end


    

    




if n_count==fig_count
    

if (plot_on==1)

% cmax=max(max_along_center);
% 

set(h_41,'Zdata',Ez')
drawnow
% 
% %figure(4)
% 
% set(fig_4_1,'YData',F_system_y(1:n))
% set(fig_4_1,'XData',t(1:n))
% set(fig_4_2,'YData',Fy_system_avg(1:n))
% set(fig_4_2,'XData',t(1:n))
% drawnow



%quiver(x(i_sample),y(j_sample),fx_avg(i_sample,j_sample)',fy_avg(i_sample,j_sample)');
% surf(x*1e6,y*1e6,fy_avg')
% % set(fig_3,'Colormap',[0.0431372560560703 0.517647087574005 0.780392169952393;0.0400560237467289 0.552100896835327 0.796078443527222;0.0369747914373875 0.586554646492004 0.811764717102051;0.033893559128046 0.621008396148682 0.82745099067688;0.0308123249560595 0.655462205410004 0.843137264251709;0.027731092646718 0.689916014671326 0.858823537826538;0.0246498603373766 0.724369764328003 0.874509811401367;0.0215686280280352 0.75882351398468 0.890196084976196;0.0184873957186937 0.793277323246002 0.905882358551025;0.0154061624780297 0.827731132507324 0.921568632125854;0.0123249301686883 0.862184882164001 0.937254905700684;0.00924369785934687 0.896638631820679 0.952941179275513;0.00616246508434415 0.931092441082001 0.968627452850342;0.00308123254217207 0.965546250343323 0.984313726425171;0 1 1;0.0555555559694767 1 1;0.111111111938953 1 1;0.16666667163372 1 1;0.222222223877907 1 1;0.277777791023254 1 1;0.333333343267441 1 1;0.388888895511627 1 1;0.444444447755814 1 1;0.5 1 1;0.555555582046509 1 1;0.611111104488373 1 1;0.666666686534882 1 1;0.722222208976746 1 1;0.777777791023254 1 1;0.833333313465118 1 1;0.888888895511627 1 1;0.944444417953491 1 1;1 1 1;1 0.923076927661896 0.923076927661896;1 0.846153855323792 0.846153855323792;1 0.769230782985687 0.769230782985687;1 0.692307710647583 0.692307710647583;1 0.615384638309479 0.615384638309479;1 0.538461565971375 0.538461565971375;1 0.461538463830948 0.461538463830948;1 0.384615391492844 0.384615391492844;1 0.307692319154739 0.307692319154739;1 0.230769231915474 0.230769231915474;1 0.15384615957737 0.15384615957737;1 0.0769230797886848 0.0769230797886848;1 0 0;0.972222208976746 0 0;0.944444417953491 0 0;0.916666686534882 0 0;0.888888895511627 0 0;0.861111104488373 0 0;0.833333313465118 0 0;0.805555582046509 0 0;0.777777791023254 0 0;0.75 0 0;0.722222208976746 0 0;0.694444417953491 0 0;0.666666686534882 0 0;0.638888895511627 0 0;0.611111104488373 0 0;0.583333313465118 0 0;0.555555582046509 0 0;0.527777791023254 0 0;0.5 0 0]);
% 
%        hold on
%        contour(x*1e6,y*1e6,er'-er_1,1,'color','red')
%        hold off
%        xlabel('x axis [{\mu}m]')
%        ylabel('y axis [{\mu}m]')
% %      c_lim=    [-2*a/eta_o 2*a/eta_o];  
% %        caxis(c_lim)
%        shading flat
% % Draw PML limes
%         line([NPML_x*dx  NPML_x*dx]*1e6, [(NPML_y)*dy (Ny-NPML_y)*dy]*1e6,'color','black');  
%         line([(Nx-NPML_x)*dx (Nx-NPML_x)*dx]*1e6, [(NPML_y)*dy (Ny-NPML_y)*dy]*1e6,'color','black');  
%         line([NPML_x*dx (Nx-NPML_x)*dx]*1e6, [(NPML_y)*dy (NPML_y)*dy]*1e6,'color','black');  
%         line([NPML_x*dx (Nx-NPML_x)*dx]*1e6, [(Ny-NPML_y)*dy (Ny-NPML_y)*dy]*1e6,'color','black'); %         
% % % Draw slab and pulse x- center of mass
% % if((x_bar_pulse(1,n)<(Nx-NPML_x)*dx)&&(y_bar_pulse(1,n)<(Ny-NPML_y)))

end


    
end

if n_count==fig_count
    n_count=0;
    
end
% 
% t(n);
n_count=n_count+1;
% x_test=x_bar_pulse(1,n)*1e6


%disp(strcat(' FDTD is_',num2str(n/Nt*100,3),' % done '))
    end
    
    
    %% save tensor data
    
  
% 
%     fx_avg(1:NPML_x+3,:)=0;
%     fy_avg(1:NPML_x+3,:)=0;
%     
%      fx_avg(Nx-NPML_x-3:Nx,:)=0;
%     fx_avg(Nx-NPML_x-3:Nx,:)=0;
%     
%        fx_avg(:,1:NPML_y+3)=0;
%     fy_avg(:,1:NPML_y+3)=0;
%     
%     fx_avg(:,Ny-NPML_y-3:Ny)=0;
%     fy_avg(:,Ny-NPML_y-3:Ny)=0;
    

    
    
fy_avg(i,j)=sum(fy_track(i,j,:),3)./N_avg;
fx_avg(i,j)=sum(fx_track(i,j,:),3)./N_avg; 


plot_mat=zeros(Nx,Ny);
plot_mat(nx1+2:nx2-2,ny1+5:ny2-5)=1;
fx_avg=fx_avg.*plot_mat;
fy_avg=fy_avg.*plot_mat;




x_norm_1=x(nx1);
x_norm_2=x(nx2);
y_norm_1=y(ny1);
y_norm_2=y(ny2);

% x and y axis for x-forc locations
%% Fx locations
% now normalize

% x and y axis for x-forc locations
%% Fx locations
% now normalize

x_fx=x(nx1-1:nx2+1)-dx/2;
x_fx=x_fx-x_norm_1;
x_fx=x_fx./(x_norm_2-x_norm_1);
% y axis for fx locations

y_fx=y(ny1-1:ny2+1);
y_fx=y_fx-y_norm_1;
y_fx=y_fx./(y_norm_2-y_norm_1);
[Xfx,Yfx]=meshgrid(x_fx,y_fx);

% solve for fluid fx locations



x_fx_cfd=([0:1:Nx_CFD-1]+.5)./(Nx_CFD-1);
y_fx_cfd=[0:1:Ny_CFD-1]./(Ny_CFD-1);
os=round((1/(Nx_CFD-1))/(1/(Nx-1)));

[XC YC]=meshgrid(x_fx_cfd,y_fx_cfd);
fx_EM=interp2(Xfx,Yfx,fx_avg(nx1-1+os:nx2+1+os,ny1-1:ny2+1),XC,YC);


%% Fy locations
% now normalize
x_fy=x(nx1-1:nx2+1);
x_fy=x_fy-x_norm_1;
x_fy=x_fy./(x_norm_2-x_norm_1);
% y axis for fx locations
y_fy=y(ny1-1:ny2+1)-dy/2;
y_fy=y_fy-y_norm_1;
y_fy=y_fy./(y_norm_2-y_norm_1);
[Xfy,Yfy]=meshgrid(x_fy,y_fy);

x_fy_cfd=([0:1:Nx_CFD-1])./(Nx_CFD-1);
y_fy_cfd=([0:1:Ny_CFD-1]+.5)./(Ny_CFD-1);

[XC YC]=meshgrid(x_fy_cfd,y_fy_cfd);
fy_EM=interp2(Xfy,Yfy,fy_avg(nx1-1:nx2+1,ny1-1:ny2+1),XC,YC);

Ez_EM=Ez(nx1:nx2,ny1:ny2);



  save(strcat('./output_data/',model,'_',source_direction,'_TE_tensor_data'),'Hy','Hx','Ez','t1x','t2x','t3x','t4x','t5x',...
        't1y','t2y','t3y','t4y','t5y','G_x','G_x_n_prev','G_y','G_y_n_prev','dx','dy','dt','Sy_avg','er','fx_avg','fy_avg')
%     

end


