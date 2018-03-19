function [ p ] = pressure_iterate_mat( u,v,dx,dy,dt,rho,p,iterations,tol ,beta)
%UNTITLED4 Summary of this function goes here
%   Detailed explanmation goes here


[Nxp,Nyp]=size(p);


%beta=1.2;

% Pressre boundary conditions with rho

rt=rho;
 rt(1,:)=1000;
 rt(end,:)=1000;
 rt(:,1)=1000;
 rt(:,end)=1000;

Q1=zeros(Nxp,Nyp);
vel_term=zeros(Nxp,Nyp);
q=2:Nxp-1;
k=2:Nyp-1;
%% PRESSURE ITERATION

% First solve for terms based on previous value before updated
% pressure


rho_xp(q,k)=rt(q,k)+rt(q+1,k);
rho_xm(q,k)=rt(q,k)+rt(q-1,k);

rho_yp(q,k)=rt(q,k)+rt(q,k+1);
rho_ym(q,k)=rt(q,k)+rt(q,k-1);
 Q1(q,k)=(1.0./(dx^2)).*(1./rho_xp(q,k)+1./rho_xm(q,k))+(1./(dy^2)).*(1./rho_yp(q,k)+1./rho_ym(q,k));
    
  
      vel_term(q,k)=(1/(2*dt)).*( (u(q,k)-u(q-1,k))/dx+(v(q,k)-v(q,k-1))/dy);
      

for i=1:iterations

 old_p=p;
for q=2:Nxp-1
    for k=2:Nyp-1
        
%         rho_xp=(rt(q,k)+rt(q+1,k));
%         rho_xm=(rt(q,k)+rt(q-1,k));
%         
%         rho_yp=(rt(q,k)+rt(q,k+1));
%         rho_ym=(rt(q,k)+rt(q,k-1));
        
        p_x_term=(1/(dx^2))*(p(q+1,k)/rho_xp(q,k)+p(q-1,k)/rho_xm(q,k));
        p_y_term=(1/(dy^2))*(p(q,k+1)/rho_yp(q,k)+p(q,k-1)/rho_ym(q,k));
      
        p(q,k)=beta.*(1/Q1(q,k))*(p_x_term+p_y_term-vel_term(q,k))+(1.0-beta)*p(q,k);
      
        
    end
    
    
end

if(max(max(abs(old_p-p))) <tol)
    break
end




end

end

