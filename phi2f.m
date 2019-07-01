function[f]=phi2f(phi)
% function[f]=phi2f(phi)
% phi2f give sthe inertial frequency in cycles per day for the latitude phi

% this is not exact because a sideral day is not exactly 24h
% f=2*2*pi/(24*60*60)*sin(phi*pi/180)*60*60*24/(2*pi);

omega = 7.2921159e-5;

f = 86400*2*omega*sin(phi*pi/180)/(2*pi);