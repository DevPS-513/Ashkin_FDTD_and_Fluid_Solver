close all
clear all
clc

addpath EM_functions
addpath material_data
addpath output_data

% FIGURE INITIALIZATION

s_d=get(0,'ScreenSize');    % Screen [0 0 width height]
sw=s_d(3);                  % are( Screen width
sh=s_d(4);                  % Screen height
fig_default_position=get(figure,'position');
fdp=fig_default_position;
fdp(4)=fdp(4)+125;
fdp(3)=fdp(3)+125;
fdp(2)=fdp(2)-125;
% figure fonts
axis_fontsize=12;
legend_fontsize=12;
Plw=1.5;                    % plotting linewidth for pulse,slab and total

set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultTextFontName', 'Arial')


% DEFINE SI UNITS

meters=1;
nm=meters*1e-9;
micrometers=1e-6;
fs=1e-15;
mu_o=4*pi*10^-7;
c=299792458;
eps_o=(1/(c*c*mu_o));
eta_o=sqrt(mu_o/eps_o);
Kg=1;
Agstr=1e-10;

% Figure Positions
% 35 pixels is the toolbar, 
% this will be different for every pc 

p1 = 	[5 sh/2 sw/3 sh/2];
p2=		[5+p1(3) sh/2 sw/3 sh/2];
p3=		[5 35 sw/3 sh/1.1];
p4=		[5+p3(3) 35 sw/3 sh/2];
p5=		[5+p4(3)+p3(3) 35 sw/3 sh/2];


p6=fdp;
p7= fdp;
p8=fdp;
p9=fdp;
p10=fdp;
p11=fdp;
p12=fdp;
p13=fdp;



sim_casev='Lossless_casner';
index_cases=[ 'n_160'];
models_fig=[ 'AB ';'MN ';'Chu';'EL ';'AMP'];
model_color=['blue  ';'red   ';'cyan  ';'black ';'yellow';'m     '];
models_fig = cellstr(models_fig);
index_cases = cellstr(index_cases);



F_folder_name='F_data';
mkdir(F_folder_name);

M_folder_name='material_data';
mkdir(F_folder_name);


force_data_n=zeros(length(index_cases),4);
force_data_n_not_norm=zeros(length(index_cases),4);
power_data_n=zeros(length(index_cases),2);
legend_data=zeros(size(index_cases));



for model_j=1:length(models_fig)
for index_j=1:length(index_cases)

    model=char(models_fig(model_j));
    index_case_n=char(index_cases(index_j));
    index_cases
    load(strcat(model,'_Workspace_data_',sim_casev,'_',index_case_n)); 
    index_cases
    

for i_n=1:Nt
    
      
    figure(1)

   surf(x(nx1:nx2),y,f_interface_x(:,:,i_n)')
   hold on
   contour(x(nx1:nx2),y,er(nx1:nx2,:)',1)
    shading flat
    view([0 90])
   % zlim([min(min(g_mech_x(nx1:nx2,:))) max(max(g_mech_x(nx1:nx2,:)))])
%     colorbar east
    grid off
    xlim([x(nx1-10) x(nx2+10)])
     title(model)
      
end
   

      
    

    
        
end
end

        
        
write_dat(strcat('./',M_folder_name,'/','Source.dat'),c./(G(:,1)*nm),abs(G(:,2)));
write_dat(strcat('./',M_folder_name,'/','Re_er.dat'),LAMBDA_drude(:)*1e9,real(er_w(1,:)));
write_dat(strcat('./',M_folder_name,'/','Im_er.dat'),LAMBDA_drude(:)*1e9,imag(er_w(1,:)));
write_dat(strcat('./',M_folder_name,'/','Re_mr.dat'),LAMBDA_drude(:)*1e9,real(mr_w(1,:)));
write_dat(strcat('./',M_folder_name,'/','Im_mr.dat'),LAMBDA_drude(:)*1e9,imag(mr_w(1,:)));
write_dat(strcat('./',M_folder_name,'/','Re_eta.dat'),LAMBDA_drude(:)*1e9,real(eta_w(1,:)));
write_dat(strcat('./',M_folder_name,'/','Im_eta.dat'),LAMBDA_drude(:)*1e9,imag(eta_w(1,:)));
        




clear all
% close all
clc
disp(' Done writing data')





