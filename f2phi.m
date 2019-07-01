function[phi]=f2phi(f)
% function[phi] = f2phi(f);
% f2phi gives the latitude corresponding to the inertial
% frequency f in cycle per day
% this inverses phi2f

f = f/(60*60*24); %inertial frequency in cycle per seconds
% omega = 2*pi/(24*60*60); % omega, rotation rate of the Earth
omega = 7.2921159e-5;

% f=inline('2*7.2722e-05*sin(x*pi/180)','x');

phi = (180/pi)*asin(2*pi*f/(2*omega));

