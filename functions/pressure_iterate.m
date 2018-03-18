function [ p ] = pressure_iterate( u,v,dx,dy,dt,rho,p,iterations,tol ,beta)
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

%% PRESSURE ITERATION

% First solve for terms based on previous value before updated
% pressure

for q=2:Nxp-1
    for k=2:Nyp-1
        
        
        
        % define density terms factor of 2 incorpateed laer.
        rho_xp=(rt(q,k)+rt(q+1,k));
        rho_xm=(rt(q,k)+rt(q-1,k));
        
        rho_yp=(rt(q,k)+rt(q,k+1));
        rho_ym=(rt(q,k)+rt(q,k-1));
        
        
      Q1(q,k)=(1.0/(dx^2))*(1/rho_xp+1/rho_xm)+(1/(dy^2))*(1/rho_yp+1/rho_ym);
        
     
      vel_term(q,k)=(1/(2*dt)).*( (u(q,k)-u(q-1,k))/dx+(v(q,k)-v(q,k-1))/dy);
      
     
        
    end
    
    
end


for i=1:iterations

 old_p=p;
for q=2:Nxp-1
    for k=2:Nyp-1
        
          rho_xp=(rt(q,k)+rt(q+1,k));
        rho_xm=(rt(q,k)+rt(q-1,k));
        
        rho_yp=(rt(q,k)+rt(q,k+1));
        rho_ym=(rt(q,k)+rt(q,k-1));
        
    p_x_term=(1/(dx^2))*(p(q+1,k)/rho_xp+p(q-1,k)/rho_xm);
      p_y_term=(1/(dy^2))*(p(q,k+1)/rho_yp+p(q,k-1)/rho_ym);
      
      p(q,k)=beta.*(1/Q1(q,k))*(p_x_term+p_y_term-vel_term(q,k))+(1.0-beta)*p(q,k);
      
        
    end
    
    
end

if(max(max(abs(old_p-p))) <tol)
    break
end




end

end

