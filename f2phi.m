function[phi]=f2phi(f);
% function[phi]=f2phi(f);
% f2phi gives the latitude corresponding to the inertial
% frequency f in cycle perday

f=f/(60*60*24); %inertial frequency in cycle per seconds
omega=2*pi/(24*60*60); % omega, rotation rate of the Earth

% f=inline('2*7.2722e-05*sin(x*pi/180)','x');

phi=asin(2*pi*f/(2*omega));
phi=phi*180/pi;

