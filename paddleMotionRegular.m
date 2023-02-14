function [chi, xi, zeta, dt, nt, s] = ...
    paddleMotionRegular(d,H,T,nT,nTi,Tdt,g,theta)
% calculates piston-type paddle motion according to a linear transfer function with raised cosine ramp applied to the first few wave periods
%
% Input data:
% d - water depth (m)
% H - wave height (m)
% T - wave period (s)
% nT - number of wave periods
% nTi - number of wave periods for start up
% Tdt - number of time steps per wave period
%
% Output data:
% chi - paddle displacement signal (m)
% xi - paddle velocity signal (m/s)
% zeta - paddle acceleration signal (m^2/s)
% dt - time increment (s)
% nt - number of time steps
%
% Author: Maciej Paprota
% Reference: M. Paprota. 2023. A twin wavemaker model for liquid sloshing in a rectangular tank. Ocean Engineering, 272, 113919

a = H/2; % wave amplitude (m)
omega = 2*pi/T; % angular wave frequency (1/s)
kappa = dispersionRelation(d,omega,g); % wave number (1/m)
s = a/(cosh(2*kappa*d)-1)*...
    (sinh(2*kappa*d)+2*kappa*d)/2; % displacement amplitude (m)
% time increment and time vector components
dt = T/Tdt; nt = nT*Tdt+1; t1 = (0:dt:nTi*T-dt)'; t2 = (nTi*T:dt:T*nT)';
% wavemaker paddle displacement function \chi(t)
chi = [s*sin(omega*t1+theta).*(-cos(omega/nTi/2*t1+theta)+1)/2; s*sin(omega*t2+theta)];
% wavemaker paddle velocity function \xi(t) = d\chi(t)/dt
xi = [s*omega*(cos(omega*t1+theta).*(1-cos(omega/nTi/2*t1))/2+ ... 
    sin(omega*t1+theta).*sin(omega/nTi/2*t1)/4/nTi); ...
    s*omega*cos(omega*t2+theta)];
% wavemaker paddle acceleration function \zeta(t) = d\xi(t)/dt
zeta = [s*omega^2/8*(sin(omega*t1+theta).*((4+1/nTi^2)*...
    cos(omega*t1/nTi/2)-4)+4/nTi*cos(omega*t1+theta).*sin(omega*t1/2/nTi)); ...
    -s*omega^2*sin(omega*t2+theta)];
end
