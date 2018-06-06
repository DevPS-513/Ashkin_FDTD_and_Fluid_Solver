

clear all
close all
clc

% Model Parameters
    c=1;
    m=1;
    k=1;

    dt_o=.01;

    sim_on=1;

% Initalize Model
    A=[0        1;
        -w_n^2 -2*zeta*w_n];

    n=1;

for n=2:100
 
    
% Initialize dt    
    dt=dt_o;
    q=q(:,n-1);

% Solve for coefficients
    K1=A*q;
    K2=A*(q+0.5*dt*A*q);
    K3=A*(q+0.5*dt*(A*( q+0.5*dt*q));
    K4=A*(q+dt*(A*+dt/2*(A*(q+dt/2*A*q))));

% Solve for the next time step
    qp=q+dt*(1/6)*[K1+2*K2+2*K3+K4];
    q=qp;    
end