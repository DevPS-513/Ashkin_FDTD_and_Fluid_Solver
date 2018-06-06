


function [fx_EM,fy_EM,eta_EM,Hz_EM,Sy_avg,fx_avg,fy_avg,x_fx_cfd,y_fx_cfd,x_fy_cfd,y_fy_cfd]=FDTD(Nx_fx,Ny_fx,Nx_fy,Ny_fy,x_CFD,y_CFD,xf_locs,yf_locs,er_low,er_high,LAMBDA,dx,beam_width,x_size,y_size,a,sim_cycles,model,plot_on,source_direction)


% Define SI units
    meters=1;
    nm=meters*1e-9;
    femptoseconds=1e-15;
    mu_o=4*pi*10^-7;
    c=299792458;
    eps_o=(1/(c*c*mu_o));
    eta_o=sqrt(mu_o/eps_o);
    Kg=1;


% FIGURE INITIALIZATION

    fig_count=150;              % time steps per plot
    n_count=0;


%% Simulation Parameters

    f=c./LAMBDA;                    % Center frequency of pulse
    T=1./f;                         % Center period

    dy=dx;
    dt=(1/2.0)*1/(c.*sqrt(1/dx^2+1/dy^2));

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
    source1_y=Ny-NPML_y-2;            % y-location of source
    a_o=1;
    S_yloc=source1_y-2;
    end

    if strcmp(source_direction,'up')
        source1_y=NPML_y+2;
        a_o=-1;
        S_yloc=source1_y+2;    
    end



%% convert CFD front percentages to EM front percentages, normalized to
% PML
nx1=NPML_x;
nx2=Nx-NPML_x+1;
ny1=NPML_y;
ny2=Ny-NPML_y+1;

% create axis on EM grid 
x_norm=x(nx1:nx2)-x(nx1);
x_norm=x_norm./(x(nx2)-x(nx1));

yf_EM_norm=interp1(xf_locs,yf_locs,x_norm,'spline');


% spline interpolate EM front



xf_EM=x(nx1)+x_norm.*(x(nx2)-x(nx1));
yf_EM=y(ny1)+yf_EM_norm.*(y(ny2)-y(ny1))-2*dy;

er_mat_1=zeros(Nx,Ny);
Nf=length(xf_EM);
for q=1:Nf
  ie=round(xf_EM(q)/dx)+1;
 je=round(yf_EM(q)/dy)+1;
 er_mat_1(ie,1:je)=1;   

 
end
 er_mat=er_mat_1;
 er=er_high.*(er_mat==1)+er_low.*(er_mat==0);
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
    n_os=5;

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

eta_EM=sqrt(er(nx1:nx2,ny1:ny2)./mr(nx1:nx2,ny1:ny2))';         % Store impedance to give to fluid simulator for plotting

    
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
sig_x=beam_width/1;
sig_x=sig_x/2;


N_x_gauss=round(gauss_x_width/dx);
x_gauss=[1:Nx]*dx;

gauss_x=gauss_create(gauss_x_avg,sig_x,x_gauss);
gauss_x=gauss_x./max(gauss_x);

gauss_x_start=gauss_x_avg-1.0*sig_x;
gauss_x_end=gauss_x_avg+1.0*sig_x;

nxg1=round(gauss_x_start/dx);
nxg2=round(gauss_x_end/dx);


    
%% Ouput Fields
% Initialize Fields and Current matricies
    Hz=zeros(Nx,Ny);
    Bz=zeros(Nx,Ny);
    Sx=zeros(Nx,Ny);
    Sy=zeros(Nx,Ny);
    Sx_flux=zeros(1,Nt);
    Sx_avg=zeros(1,Nt);
    Hz_at_x=zeros(Nx,Ny);
    Hz_at_y=zeros(Nx,Ny);    
    Bz_at_x=zeros(Nx,Ny);
    Bz_at_y=zeros(Nx,Ny);   
    
    Ex=zeros(Nx,Ny);
    Ey=zeros(Nx,Ny);
    Dx=zeros(Nx,Ny);
    Dy=zeros(Nx,Ny);

    % Dispersive current densities
    Jxd=zeros(Nx,Ny);
    Jyd=zeros(Nx,Ny);
    Jmzd=zeros(Nx,Ny);
        
    % Static current densities
    Jxs=zeros(Nx,Ny);
    Jys=zeros(Nx,Ny);
    Jmzs=zeros(Nx,Ny);
    
    % Dispersive polarization and magnetization.
    Pxd=zeros(Nx,Ny);
    Pyd=zeros(Nx,Ny);
    Mzd=zeros(Nx,Ny);
    
    % Static Polarization and magnetizaton
    Pxs=zeros(Nx,Ny);
    Pys=zeros(Nx,Ny);
    Mzs=zeros(Nx,Ny);
    
    % Total polarization and currents
    
    Px=zeros(Nx,Ny);
    Py=zeros(Nx,Ny);
    Mz=zeros(Nx,Ny);
    Jmz=zeros(Nx,Ny);
    Jex=zeros(Nx,Ny);
    Jey=zeros(Nx,Ny);
    
%% Output Variables

 Sy_avg=zeros(1,Nt);
    

% Total variables
    Tx=zeros(Nx,Ny);                                % Tx component
    Ty=zeros(Nx,Ny);                                % Ty component

    g_mech_x=zeros(Nx,Ny);                         %Momentum density
    g_mech_y=zeros(Nx,Ny);                         %Momentum density


% Pulse_Output_Variables
    W=zeros(Nx,Ny);                     %[J/m^3]        Energy Density

    % pulse x-variables
    G_x=zeros(Nx,Ny);                   %[Kg/(m^2*s)]   Momentum Density
    G_x_n_prev=0;

    % y-variables
    G_y=zeros(Nx,Ny);                   %[Kg/(m^2*s)]   Momentum Density
    G_y_n_prev=0;

       
    fx=zeros(Nx,Ny);
    fy=zeros(Nx,Ny);


    N_avg=round(5*T/dt);
    avg_count=1;
    fx_avg=zeros(Nx,Ny);
    fy_avg=zeros(Nx,Ny);

    fx_track=zeros(Nx,Ny,N_avg);
    fy_track=zeros(Nx,Ny,N_avg);


 

%% Update Coefficients
% Ex UPDATE COEFFCIENTS

C1=eps_o*(er(2:end,2:end)-1);

P1=(1/dt-gamma_e/2)./(1/dt+gamma_e/2);      
P2=eps_o*omega_e.^2./(2.*(1/dt+gamma_e/2));
den_x=(sigma_e/2+eps_o*er/dt+1/2.*P2);
X1=(eps_o*er/dt+sigma_e/2+P2/2);
X2=(eps_o*er/dt-sigma_e/2-P2/2);

X3=X2./X1;
X4=0.5.*((P1+1)./X1);
X5=(1./(dy.*X1));

P1=P1(2:end,2:end);
P2=P2(2:end,2:end);

X3=X3(2:end,2:end);
X4=X4(2:end,2:end);
X5=X5(2:end,2:end);

% Hz Update Coefficients

C2=mu_o*(mr-1);

M1=(1/dt-gamma_m./2)./(gamma_m./2+1/dt);                                    
M2=mu_o.*omega_m.^2./(2.*(gamma_m./2+1/dt));                                   

A1=(mu_o*mr./dt+sigma_m./2+1/2.*M2);
A2=(mu_o*mr./dt-sigma_m./2-1/2.*M2);

A3=A2./A1;

A4=(1./(2.*A1));
M3=M1+1;
A5=1./A1;

M1=M1(1:end-1,1:end-1);
M2=M2(1:end-1,1:end-1);
M3=M3(1:end-1,1:end-1);
A3=A3(1:end-1,1:end-1);
A4=A4(1:end-1,1:end-1);
A5=A5(1:end-1,1:end-1);


% Ey UPDATE COEFFICIENTS

P1Y=(1/dt-gamma_e/2)./(1/dt+gamma_e/2);
P2Y=eps_o*omega_e.^2./(2*(1/dt+gamma_e/2));
den_y=(sigma_e/2+eps_o*er/dt+1/2.*P2Y);

Y1=(eps_o*er/dt-sigma_e/2-P2Y/2)./den_y;
Y2=(P1Y+1)./(2*(den_y));
Y3=((1/dy)*(1./(den_y)));

P1Y=P1Y(2:end,2:end);
P2Y=P2Y(2:end,2:end);

Y1=Y1(2:end,2:end);
Y2=Y2(2:end,2:end);
Y3=Y3(2:end,2:end);

%% Indexes

i=4:Nx-4;
j=4:Ny-4;
i_W=i;
j_W=[4:1:source1_y-4, source1_y+4:1:Ny-4]; 


if plot_on==1
    figure(2)
    surf(([nx1:nx2]-nx1)./(nx2-nx1),([ny1:ny2]-ny1)./(ny2-ny1),er(nx1:nx2,ny1:ny2)')
    hold on
    % plot3(xf_locs,yf_locs,10.*ones(size(xf_locs)),'-r')
    % hold on
    plot3(x_norm,yf_EM_norm,10.*ones(size(x_norm)),'g')
    xlim([-.2 1.2 ])
    ylim([-.2 1.2])
    view([ 0 90])
    shading flat
    title(' index n in FDTD')

    hold off
    figure(4)
    h_41=surf(x,y,Hz');
    shading flat
    view([0 90])
    hold off
    figure(5)
    surf(x*1e6,y*1e6,sqrt(mr./er)')
    shading flat
    view([0 90])
    title(' \eta in FDTD')
    hold on
    plot3(xf_EM*1e6,yf_EM*1E6,10.*a^2.*ones(size(xf_EM)),'linewidth',1.8,'color','g')
    line([NPML_x*dx  NPML_x*dx]*1e6, [(NPML_y)*dy (Ny-NPML_y)*dy]*1e6,[10 10],'color','black');  
    line([(Nx-NPML_x)*dx (Nx-NPML_x)*dx]*1e6, [(NPML_y)*dy (Ny-NPML_y)*dy]*1e6,[10 10],'color','black');  
    line([NPML_x*dx (Nx-NPML_x)*dx]*1e6, [(NPML_y)*dy (NPML_y)*dy]*1e6,[10 10],'color','black');  
    line([NPML_x*dx (Nx-NPML_x)*dx]*1e6, [(Ny-NPML_y)*dy (Ny-NPML_y)*dy]*1e6,[10 10],'color','black'); % 
    hold off


end

% 



for n= 1:Nt
      

% Store Previous values

    Ex_n_prev=Ex(2:Nx,2:Ny);
    Dx_n_prev=Dx(2:Nx,2:Ny);   
    Jxd_n_prev=Jxd(2:Nx,2:Ny);      
    Jxs_n_prev=Jxs(2:Nx,2:Ny);      

    Ey_n_prev=Ey(2:Nx,2:Ny);
    Dy_n_prev=Ey(2:Nx,2:Ny);
    Jyd_n_prev=Jyd(2:end,2:end);   
    Jys_n_prev=Jys(2:end,2:end);   


    Hz_n_prev=Hz;
    Bz_n_prev=Bz;
    Jmz_n_prev_fr=Jmz;
    Jmzd_n_prev=Jmzd;    
              

%% Update Hz (n+1/2)

    E_term=-1*(1/dx).*(Ey(2:Nx,1:Ny-1)-Ey(1:Nx-1,1:Ny-1))...
            +(1/dy).*(Ex(1:Nx-1,2:Ny)-Ex(1:Nx-1,1:Ny-1));

    Hz(1:Nx-1,1:Ny-1)=  ((A3).*Hz(1:Nx-1,1:Ny-1)...    
                        -(A4).*((M3).*Jmzd(1:Nx-1,1:Ny-1))...
                        +(A5).*E_term);

    Jmzd(1:Nx-1,1:Ny-1)=M1.*Jmzd(1:Nx-1,1:Ny-1)+M2.*(Hz(1:Nx-1,1:Ny-1)...
                        +Hz_n_prev(1:Nx-1,1:Ny-1));            

    Jmzs(1:Nx-1,1:Ny-1)=-1*Jmzs(1:Nx-1,1:Ny-1)+(2/dt).*C2(1:Nx-1,1:Ny-1).*(Hz(1:Nx-1,1:Ny-1)-Hz_n_prev(1:Nx-1,1:Ny-1));
    Jmz(1:Nx-1,1:Ny-1)=Jmzs(1:Nx-1,1:Ny-1)+Jmzd(1:Nx-1,1:Ny-1);

        
    Mzd(1:Nx-1,1:Ny-1)=Mzd(1:Nx-1,1:Ny-1)+(dt/2).*(Jmzd(1:Nx-1,1:Ny-1)+Jmzd_n_prev(1:Nx-1,1:Ny-1));
    Mzs(1:Nx-1,1:Ny-1)=C2(1:Nx-1,1:Ny-1).*Hz(1:Nx-1,1:Ny-1);
    Mz=Mzd+Mzs;
    Bz=mu_o.*Hz+Mz;
             
        
        eta_r=sqrt(mr(round(Nx/2),source1_y)./er(round(Nx/2),source1_y));                                                    % Scale Hz amplitude for one way source
        Hz(nxg1:nxg2,source1_y)=a_o*gauss_x(nxg1:nxg2)'.*(a/(2*eta_o*eta_r)).*sin(2*pi*f*(t(n)))+ Hz(nxg1:nxg2,source1_y);            % X- Guassian %         
    


if (strcmp(model,'MN'))
        %[W]=calculate_W_MN(c^2,i_W,j_W,Dx,Dx_n_prev,Dy,Dy_n_prev,Bz,Bz_n_prev,dx,dy,dt,W);
        % turn off W calculation for speed
        %    W(source1_x-1:source1_x+1,gi)=0;                                   % Remove energy at source location

        Bz_at_y(i,j)=(1/4)*(Bz(i,j)+Bz(i-1,j)+Bz_n_prev(i,j)+Bz_n_prev(i-1,j));
        Bz_at_x(i,j)=(1/4)*(Bz(i,j)+Bz(i,j-1)+Bz_n_prev(i,j)+Bz_n_prev(i,j-1));

      % Bring Bz in space to Dx, but not in time

        Bz_at_x(i,j)=(1/2)*(Bz(i,j)+Bz(i,j-1));


        G_x_n_prev=G_x;
        G_x(i,j)=Dy(i,j).*Bz_at_y(i,j);                  % G_MN=DXB
        Sx(i,j)=(Dy(i,j).*Bz_at_y(i,j)).*c^2;

        G_y_n_prev=G_y;
        G_y(i,j)=-1*Dx(i,j).*Bz_at_x(i,j);   % Ty is at Ex i think    

        % Tx for momentum continuity equation nabla*T
        [Tx ] = Calculate_Tx_MN( i,j,Ex,Ey,Dx,Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );
       % [Ty ] = Calculate_Ty_MN( i,j,Ex,Ey,Dx,Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );
        [Ty,t1,t2,t3,t4 ] = Calculate_Ty_MN( i,j,Ex,Ex_n_prev,Ey,Ey_n_prev,Dx,Dx_n_prev,...
            Dy,Dy_n_prev,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );

    
end


if (strcmp(model,'Chu'))    
        %[W]=calculate_W_AB_v2(i_W,j_W,Ex,Ex_n_prev,Ey,Ey_n_prev,Hz,Hz_n_prev,dx,dy,dt,W);  
        % W(source1_x-1:source1_x+1,gi)=0;

        Hz_at_y(i,j)=(1/4)*(Hz(i,j)+Hz(i-1,j)+Hz_n_prev(i,j)+Hz_n_prev(i-1,j)); 
        Hz_at_x(i,j)=(1/4)*(Hz(i,j)+Hz(i,j-1)+Hz_n_prev(i,j)+Hz_n_prev(i,j-1));    

        G_x_n_prev=G_x;        
        G_x(i,j)=eps_o*mu_o.*(Ey(i,j).*Hz_at_y(i,j)); % Chu, EL and AB  
        Sx(i,j)=(Ey(i,j).*Hz_at_y(i,j));

        [ Tx] = Calculate_Tx_Chu(i,j,Ex,Ey,Dx,Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );   

        G_y_n_prev=G_y;
        G_y(i,j)=-1*eps_o*mu_o.*(Ex(i,j).*Hz_at_x(i,j));               
        [Ty,t1,t2,t3,t4 ] = Calculate_Ty_Chu( i,j,Ex,Ey,Dx,Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );
end


if (strcmp(model,'EL'))    
        %  [W]=calculate_W_AB_v2(i_W,j_W,Ex,Ex_n_prev,Ey,Ey_n_prev,Hz,Hz_n_prev,dx,dy,dt,W);
        %   W(source1_x-1:source1_x+1,gi)=0;

        Hz_at_y(i,j)=(1/4)*(Hz(i,j)+Hz(i-1,j)+Hz_n_prev(i,j)+Hz_n_prev(i-1,j)); 
        Hz_at_x(i,j)=(1/4)*(Hz(i,j)+Hz(i,j-1)+Hz_n_prev(i,j)+Hz_n_prev(i,j-1));
        G_x_n_prev=G_x;        
        G_x(i,j)=eps_o*mu_o.*(Ey(i,j).*Hz_at_y(i,j)); % Chu, EL and AB    
        [Tx] = Calculate_Tx_EL(i,j,Ex,Ey,Dx,Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );   
        Sx(i,j)=(Ey(i,j).*Hz_at_y(i,j));
        
        G_y_n_prev=G_y;
        G_y(i,j)=-1*eps_o*mu_o.*(Ex(i,j).*Hz_at_x(i,j));                  % 
        % Ty is derivatives about x location 
        % Ty is at i+1/2, and in time at B
        % [Ty,t1,t2,t3,t4 ] = Calculate_Ty_EL_v3( i,j,Ex,Ex_n_prev,Ey,Ey_n_prev,Dx,Dx_n_prev,...
        % Dy,Dy_n_prev,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );  
    
        [Ty,t1,t2,t3,t4 ] = Calculate_Ty_EL( i,j,Ex,Ey,Dx,...
         Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );  
end

if (strcmp(model,'AB'))
    
       % [W]=calculate_W_AB_v2(i,j,Ex,Ex_n_prev,Ey,Ey_n_prev,Hz,Hz_n_prev,dx,dy,dt,W); 
       %  W(source1_x-1:source1_x+1,gi)=0;

        Hz_at_y(i,j)=(1/4)*(Hz(i,j)+Hz(i-1,j)+Hz_n_prev(i,j)+Hz_n_prev(i-1,j)); 
        Hz_at_x(i,j)=(1/4)*(Hz(i,j)+Hz(i,j-1)+Hz_n_prev(i,j)+Hz_n_prev(i,j-1));
        G_x_n_prev=G_x;        
        G_x(i,j)=eps_o*mu_o.*(Ey(i,j).*Hz_at_y(i,j)); % Chu, EL and AB  
        
        [ Tx] = Calculate_Tx_AB(i,j,Ex,Ey,Dx,Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );   
        Sx(i,j)=(Ey(i,j).*Hz_at_y(i,j));

        G_y_n_prev=G_y;
        G_y(i,j)=-1*eps_o*mu_o.*(Ex(i,j).*Hz_at_x(i,j));                  % 
        %[Ty,t1,t2,t3,t4 ] = Calculate_Ty_AB( i,j,Ex,Ey,Dx,Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );

         % MN and AB have the same stress tensor for all materials
         % considered
         
        [Ty,t1,t2,t3,t4 ] = Calculate_Ty_MN( i,j,Ex,Ex_n_prev,Ey,Ey_n_prev,Dx,Dx_n_prev,...
            Dy,Dy_n_prev,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );
    
end

if (strcmp(model,'AMP'))
    
       % [W]=calculate_W_MN(eps_o*c^2,i_W,j_W,Ex,Ex_n_prev,Ey,Ey_n_prev,Bz,Bz_n_prev,dx,dy,dt,W);  
     %    W(source1_x-1:source1_x+1,gi)=0;

        Bz_at_y(i,j)=(1/4)*(Bz(i,j)+Bz(i-1,j)+Bz_n_prev(i,j)+Bz_n_prev(i-1,j));   
        Bz_at_x(i,j)=(1/4)*(Bz(i,j)+Bz(i,j-1)+Bz_n_prev(i,j)+Bz_n_prev(i,j-1));
        G_x_n_prev=G_x;        
        G_x(i,j)=eps_o*(Ey(i,j).*Bz_at_y(i,j));
        [Tx] = Calculate_Tx_AMP(i,j,Ex,Ey,Dx,Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );  
        G_y_n_prev=G_y;
        G_y(i,j)=-1*eps_o*(Ex(i,j).*Bz_at_x(i,j));                 % 
        [Ty,t1,t2,t3,t4 ] = Calculate_Ty_AMP( i,j,Ex,Ey,Dx,Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );
       
        
        Sx(i,j)=eps_o*(Ey(i,j).*Bz_at_y(i,j)).*c^2;
        %Sx(i,j)=mu_o.*(Dy(i,j).*Hz_at_y(i,j)).*c^2;
        G_x(i,j)=Sx(i,j)/c^2;    

end


 
%% Update Ex,Ey (n+1)

% Update Ex
    Ex(2:Nx,2:Ny)=(X3).*Ex(2:Nx,2:Ny)-X4.*Jxd(2:Nx,2:Ny)+X5.*(Hz(2:Nx,2:Ny)-Hz(2:Nx,1:Ny-1));
    
    Jxd(2:Nx,2:Ny)=P1.*Jxd(2:Nx,2:Ny)+P2.*(Ex(2:Nx,2:Ny)+Ex_n_prev);
    Pxd(2:Nx,2:Ny)=Pxd(2:Nx,2:Ny)+(dt/2).*(Jxd(2:Nx,2:Ny)+Jxd_n_prev);

    Jxs(2:Nx,2:Ny)=-Jxs(2:Nx,2:Ny)+(2/dt)*C1.*(Ex(2:Nx,2:Ny)-Ex_n_prev);
    Pxs(2:Nx,2:Ny)=Pxs(2:Nx,2:Ny)+(dt/2)*(Jxs(2:Nx,2:Ny)+Jxs_n_prev);

    Jex(2:Nx,2:Ny)=Jxs(2:Nx,2:Ny)+Jxd(2:Nx,2:Ny);    
    Px(2:Nx,2:Ny)=Pxd(2:Nx,2:Ny)+Pxs(2:Nx,2:Ny);
    Dx=eps_o.*Ex+Px;

% UPDATE Ey


    Ey(2:Nx,2:Ny)=Y1.*Ey(2:Nx,2:Ny)-Y2.*Jyd(2:Nx,2:Ny)-Y3.*(Hz(2:Nx,2:Ny)-Hz(1:Nx-1,2:Ny));

    Jyd(2:Nx,2:Ny)=P1Y.*Jyd(2:Nx,2:Ny)+P2Y.*(Ey(2:Nx,2:Ny)+Ey_n_prev);
    Pyd(2:end,2:end)=Pyd(2:end,2:end)+(dt/2)*(Jyd(2:end,2:end)+Jyd_n_prev);

    Jys(2:Nx,2:Ny)=-Jys(2:Nx,2:Ny)+(2/dt)*C1.*(Ey(2:Nx,2:Ny)-Ey_n_prev);
    Pys(2:end,2:end)=Pys(2:end,2:end)+(dt/2)*(Jys(2:end,2:end)+Jys_n_prev);

    Jey(2:Nx,2:Ny)=Jys(2:Nx,2:Ny)+Jyd(2:Nx,2:Ny);
    Py(2:end,2:end)=Pyd(2:end,2:end)+Pys(2:end,2:end);
    Dy=eps_o.*Ey+Py;
 

    Ex(nxg1:nxg2,source1_y)=gauss_x(nxg1:nxg2)'.*(a/(2)).*sin(2*pi*f*(t(n)))+ Ex(nxg1:nxg2,source1_y);            % Ex for one way source        


%% Slab Calculations


    g_mech_x_prev=g_mech_x;
    g_mech_y_prev=g_mech_y;


    g_mech_x(i,j)=g_mech_x(i,j)-dt.*Tx(i,j)-1*(G_x(i,j)-G_x_n_prev(i,j));
    g_mech_y(i,j)=g_mech_y(i,j)-dt.*Ty(i,j)-1*(G_y(i,j)-G_y_n_prev(i,j));
    
    
    fx(i,j)=(g_mech_x(i,j)-g_mech_x_prev(i,j))/dt;
    fy(i,j)=(g_mech_y(i,j)-g_mech_y_prev(i,j))/dt;

    fx(:,source1_y-6:source1_y+6)=0;
    fy(:,source1_y-6:source1_y+6)=0;
    
     
% Start a counter, it will store the instantaneous fx and fy
% over N_avg time steps
    if avg_count<=N_avg
         fx_track(:,:,avg_count)=fx;
         fy_track(:,:,avg_count)=fy;
         avg_count=avg_count+1;
             if avg_count==N_avg
                avg_count=1; 
             end
     end



    if n_count==fig_count

    % Plot the electric field 
            if (plot_on==1)


            set(h_41,'Zdata',Hz')
            drawnow


            end



    end

    if n_count==fig_count
        n_count=0;

    end

n_count=n_count+1;

 end

    
% Average the fx and fy force densities over N_avg samples taken
% these are large 2d matricies
fy_avg(i,j)=sum(fy_track(i,j,:),3)./N_avg;
fx_avg(i,j)=sum(fx_track(i,j,:),3)./N_avg; 



% Create start and stop points for FDTD grid
% x_norm_1 will be where the first point on the x axis is considered
% The PML is ignored

x_norm_1=x(NPML_x);
x_norm_2=x(Nx-NPML_x+1);
y_norm_1=y(NPML_x);
y_norm_2=y(Ny-NPML_y+1);

% x and y axis for x-forc locations
%% Fx locations, normalize the average fx to fit onto the fluid grid

x_fx=x(nx1-1:nx2+1)-dx/2;
x_fx=x_fx-x_norm_1;
x_fx=x_fx./(x_norm_2-x_norm_1);
% y axis for fx locations

y_fx=y(ny1-1:ny2+1);
y_fx=y_fx-y_norm_1;
y_fx=y_fx./(y_norm_2-y_norm_1);
[Xfx,Yfx]=meshgrid(x_fx,y_fx);

% solve for fluid fx locations
x_fx_cfd=[0:1:Nx_fx-1]./(Nx_fx-1);
y_fx_cfd=[0:1:Ny_fx-1]./(Ny_fx-1)+round(0.5/(Ny_fx-1));

[XC YC]=meshgrid(x_fx_cfd,y_fx_cfd);
fx_EM=interp2(Xfx,Yfx,fx_avg(nx1-1:nx2+1,ny1-1:ny2+1),XC,YC);
fx_EM(isnan(fx_EM))=0;

%% Fy locations, normalize the average Fy to fit onto the fluid grid

x_fy=x(nx1-1:nx2+1);
x_fy=x_fy-x_norm_1;
x_fy=x_fy./(x_norm_2-x_norm_1);
% y axis for fx locations
y_fy=y(ny1-1:ny2+1)-dy/2;
y_fy=y_fy-y_norm_1;
y_fy=y_fy./(y_norm_2-y_norm_1);
[Xfy,Yfy]=meshgrid(x_fy,y_fy);


x_fy_cfd=[0:1:Nx_fy-1]./(Nx_fy-1)+round(0.5/(Nx_fy-1));
y_fy_cfd=[0:1:Ny_fy-1]./(Ny_fy-1);

[XC YC]=meshgrid(x_fy_cfd,y_fy_cfd);



fy_EM=interp2(Xfy,Yfy,fy_avg(nx1-1:nx2+1,ny1-1:ny2+1),XC,YC);
fy_EM(isnan(fy_EM))=0;

Hz_EM=Hz(nx1:nx2,ny1:ny2);

f_max=max([ max(max(abs(fy_EM))) max(max(abs(fx_EM))) ]);
fx_EM=fx_EM./f_max;
fy_EM=fy_EM./f_max;



end


