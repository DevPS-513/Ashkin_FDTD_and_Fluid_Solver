function [ G ] = gauss_create( mu,sig,x )
% returns a matrix of frequency values of coefficients 
% for a gaussian power spectrum centered at f_c



G=(1/sig*sqrt(2*pi)).*exp(-1/2*((x-mu)/(sig)).^2);






end

