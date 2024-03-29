function [kappa, eta_hat, phi_hat] = ...
    wavemakerAnalytical(l,d,chi_hat,sigma,t,I,g)
% analytical solution to linear mechanically-generated gravity waves in a closed wave flume of constant depth generated by a piston-type wavemaker moving according to \chi = \chi_hat*\sin(\sigma t) and reflected at a vertical wall
% 
% Input data:
% l - wave flume length (m)
% d - water depth (m)
% chi_hat - wavemaker paddle motion amplitude (m)
% sigma - wavemaker paddle motion frequency (rad/s)
% t - time domain vector (s)
% I - number of solution eigenvalues
% g - gravitational acceleration (m/s^2)
%
% Output data:
% kappa - solution eigenvalues (rad/m)
% eta_hat - free-surface amplitudes (m) 
% phi_hat - velocity potential amplitudes (m^2/s)
%
% Author: Maciej Paprota
% Reference: M. Paprota. 2023. A twin wavemaker model for liquid sloshing in a rectangular tank. Ocean Engineering, 272, 113919

kappa = (1:I)*pi/l; % solution eigenvalues (rad/m)
omega = sqrt(g*kappa.*tanh(kappa*d)); % solution frequences (rad/s)
alpha = -2*chi_hat*sigma^2*tanh(kappa*d)./kappa/l;
eta_hat = (sin(t*omega).*omega/sigma-sin(t*sigma)).*...
    alpha./(sigma^2-omega.^2); % eq. (13)
phi_hat = (cos(t*sigma)*sigma^2-cos(t*omega).*omega.^2).*...
    alpha*g./omega.^2/sigma./(omega.^2-sigma^2); % eq. (14)
eta_hat(:,omega==sigma) = -(sigma*cos(t*sigma).*t+sin(t*sigma))*...
    alpha(omega==sigma)/2/sigma^2; % resonance eq. (13)
phi_hat(:,omega==sigma) = (sigma*sin(t*sigma).*t-2*cos(sigma*t))*...
    alpha(omega==sigma)*g/2/sigma^3; % resonance eq. (14)
eta_hat = [sin(sigma*t)*chi_hat*d/l eta_hat]; % + eq. (12)
end
