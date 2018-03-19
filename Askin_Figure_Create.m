clear all
close all
clc

addpath output_data
addpath EM_functions
addpath functions
addpath material_data


s_d=get(0,'ScreenSize');  % Screen [0 0 width height]
sw=s_d(3);                % Screen width
sh=s_d(4);                % Screen height

models=cellstr(['AB ';'MN ']);
% models=cellstr(['MN'])
s_cases=cellstr(['up']);
mkdir('./output_data/')
mkdir('./Figures')



    
    for case_j=1:length(s_cases)
        
    
    source_direction=char(s_cases(case_j));
       AB=load(strcat(pwd,'./output_data/','AB','_',source_direction));
        MN=load(strcat(pwd,'./output_data/','MN','_',source_direction));
  EL=load(strcat(pwd,'./output_data/','EL','_',source_direction));
    AMP=load(strcat(pwd,'./output_data/','AMP','_',source_direction));
dx=AB.data.parameters(1);       
dy=AB.data.parameters(2);   
Nx=AB.data.parameters(3);       
Ny=AB.data.parameters(4); 
dt=AB.data.parameters(5); 
Nt=AB.data.parameters(6); 
%fig_count=data.parameters(7); 
fig_count=1;

for iq=400:Nt
 
   
    xf_AB=AB.data.(strcat('AB','_xf_',source_direction,'_',num2str(iq)));
    yf_AB=AB.data.(strcat('AB','_yf_',source_direction,'_',num2str(iq)));
      dx=AB.data.parameters(1); 
      dy=dx;
     eta_AB=AB.data.(strcat('AB','_eta_EM_',source_direction,'_',num2str(iq)));
     Hz_AB=AB.data.(strcat('AB','_Hz_EM_',source_direction,'_',num2str(iq)));
r_AB=AB.data.(strcat('AB','_r_',source_direction,'_',num2str(iq)));

fx_EM_AB=AB.data.(strcat('AB','_fx_EM_',source_direction,'_',num2str(iq)));
 fy_EM_AB=AB.data.(strcat('AB','_fy_EM_',source_direction,'_',num2str(iq)));   
 
 
    xf_MN=MN.data.(strcat('MN','_xf_',source_direction,'_',num2str(iq)));
    yf_MN=MN.data.(strcat('MN','_yf_',source_direction,'_',num2str(iq)));
    
    xf_EL=EL.data.(strcat('EL','_xf_',source_direction,'_',num2str(iq)));
    yf_EL=EL.data.(strcat('EL','_yf_',source_direction,'_',num2str(iq)));
    
    xf_AMP=AMP.data.(strcat('AMP','_xf_',source_direction,'_',num2str(iq)));
    yf_AMP=AMP.data.(strcat('AMP','_yf_',source_direction,'_',num2str(iq)));
  
    
%     fx=data.(strcat(model,'_fx_',source_direction,'_',num2str(iq*dt*1E6,'%u')));
%     fy=data.(strcat(model,'_fy_',source_direction,'_',num2str(iq*dt*1E6,'%u')));
%     fx_EM=data.(strcat(model,'_fx_EM_',source_direction,'_',num2str(iq*dt*1E6,'%u')));
%     fy_EM=data.(strcat(model,'_fy_EM_',source_direction,'_',num2str(iq*dt*1E6,'%u')));

%% FIGURE INITIALIZE
if iq==400
    f_1=figure(1);
    h_11=fill(xf_AB(3:end-2),yf_AB(3:end-2),[136/255 216/255 247/255],'linewidth',2);
    hold on
    h_12=plot(xf_MN(3:end-2),yf_MN(3:end-2),'--r','linewidth',3);
    xlim([-dx/2 Nx*dx])
    ylim([-dy/2 Ny*dy])
%     legend('AB','MN')
    set(f_1,'position',[0 sh/2-35 sw/3 sh/2])
    axis off
    axis equal
 
    
    
    f_2=figure(2);
    h_21=plot(xf_AB(2:end-1),yf_AB(2:end-1),'linewidth',2.5);
    hold on
    h_22=plot(xf_EL(2:end-1),yf_EL(2:end-1),'--black','linewidth',2);
    xlim([-dx/2 Nx*dx])
    ylim([-dy/2 Ny*dy])
    legend('AB','EL')
    set(f_2,'position',[sw/3 sh/2-35 sw/3 sh/2])
    
    f_3=figure(3);
    h_31=plot(xf_AB(2:end-1),yf_AB(2:end-1),'linewidth',2.5);
    hold on
    h_32=plot(xf_AMP(2:end-1),yf_AMP(2:end-1),'--y','linewidth',2);
    xlim([-dx/2 Nx*dx])
    ylim([-dy/2 Ny*dy])
    legend('AB','AMP')
    set(f_3,'position',[2*sw/3 sh/2-35 sw/3 sh/2])
    
else
       if iq==1000
  
        [Nx_EM,Ny_EM]=size(Hz_AB);
        dx_em=1./(Nx_EM-1);
        dy_em=1./(Ny_EM-1);
        x_em=[0:1:Nx_EM-1]*dx_em;
        y_em=[0:1:Ny_EM-1]*dy_em;
  
        [Nx Ny]=size(r_AB);
        Nx=Nx-1;
        Ny=Ny-1;
        

        max_f=max(max(sqrt(fx_EM_AB(1:Nx-1,1:Ny-1).^2+fy_EM_AB(1:Nx-1,1:Ny-1).^2)));
        fx_EM=fx_EM_AB./max_f;
        fy_EM=fy_EM_AB./max_f;
        
        fx_EM=fx_EM*20/Nx;
        fy_EM=fy_EM*20/Ny;

        fx_EM=fx_EM(:,1:end-1);
        fy_EM=fy_EM(1:end-1,:);
        
        
        % being all values to middle
        
        

        
        
        x=[0:1:Nx-1]./(Nx-1);
        y=[0:1:Ny-1]./(Ny-1);
        
        for i=1:2:Nx-1
            for j=1:2:Ny-1
                
                
                
                f_mag=sqrt(fx_EM(i,j)^2+fy_EM(i,j)^2)*(Ny/20);
                if(f_mag>.35)
                    
                    fx_EM(i,j)=.5*(fx_EM(i,j)+fx_EM(i,j+1));
                    fy_EM(i,j)=.5*(fy_EM(i,j)+fy_EM(i+1,j));
                    
                    figure(10)
                    hold on
                    h=quiver(x(i),y(j),fx_EM(i,j),fy_EM(i,j));
                    set(h,'color','r');
                                     hold on                   
                    
                    
                end
                
                
                
            end
            
        end
        
     figure(1)
        hold off
        
  %     saveas(f_1,'./Figures/Ashkin_geom','png'); 
    end
    
    
    set(h_11,'XData',xf_AB(3:end-3),'YData',yf_AB(3:end-3))
        set(h_12,'XData',xf_MN(2:end-1),'YData',yf_MN(2:end-1))
    set(h_21,'XData',xf_AB(2:end-1),'YData',yf_AB(2:end-1))
 set(h_22,'XData',xf_EL(2:end-1),'YData',yf_EL(2:end-1))
 
     set(h_31,'XData',xf_AB(2:end-1),'YData',yf_AB(2:end-1))
 set(h_32,'XData',xf_AMP(2:end-1),'YData',yf_AMP(2:end-1))
 
drawnow
 
end
if iq>250
%     pause(.1);
end
    iq
end
    
    

%  clf

end


       
        
        
        
        
        
        
        

