function [ er] = drude_calc( Gamma,wp,e_inf,Lambda )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

f=299792458/Lambda;

w=2*3.14159265358979*f;

T=Gamma;


real=e_inf-(wp^2)/(w^2+T^2);


imag=((wp^2)*T)/(w^3+w*T^2);

er=real+j*imag;
end

