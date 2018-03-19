function [ G ] = gauss_pulse( tau, f_c,FWHM,N,f_min,f_max )
% returns a matrix of frequency values of coefficients 
% for a gaussian power spectrum centered at f_c
df=(f_max-f_min)/N;
f=[f_min:df:f_max];
sig=FWHM/(2*sqrt(2*log(2)));
G=zeros(length(f),2);
G(:,2)=exp(-i*2*pi*f*tau).*exp(-1*((f-f_c)./sig).^2);
G(:,1)=f;
end

