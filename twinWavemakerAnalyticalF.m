function [lambda, phi_hat] = ...
    twinWavemakerAnalyticalF(l,d,chi_hat,sigma,t,K,g)
% analytical solution to linear sloshing in a closed tank of constant depth generated by twin piston-type wavemakers moving according to \chi = \chi_hat*\sin(\sigma t).
% Faltinsen (1978) version
%
% Input data:
% l - wave flume length (m)
% d - water depth (m)
% chi_hat - wavemaker paddle motion amplitude (m)
% sigma - wavemaker paddle motion frequency (rad/s)
% t - time domain vector (s)
% K - number of solution eigenvalues
% g - gravitational acceleration (m/s^2)
%
% Output data:
% lambda - solution eigenvalues (rad/m)
% phi_hat - velocity potential amplitudes (m^2/s)
%
% Author: Maciej Paprota
% Reference: M. Paprota. 2023. A twin wavemaker model for liquid sloshing in a rectangular tank. Ocean Engineering, 272, 113919

lambda = (2*(1:K)-1)*pi/l; % solution eigenvalues (rad/m)
omega = sqrt(g*lambda.*tanh(lambda*d)); % solution frequences (rad/s)
phi_hat = (cos(t*sigma)*sigma^2-cos(t*omega).*omega.^2).*...
    (-1).^(1:K).*4*chi_hat*sigma./(sigma^2-omega.^2)./lambda.^2/l; % amplitudes in eq. (32)
end
