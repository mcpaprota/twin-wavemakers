function eta = freeSurfaceElevation(eta_hat,kappa,x)
% calculates free-surface elevation based on amplitudes and eigenvalues of cosine expansion
%
% Input data:
% eta_hat - solution amplitudes (m)
% kappa(lambda) - solution aigenvalues (1/m)
% x - longitudinal coordinates (m)
%
% Output data:
% eta - free surface elevation (m)
%
% Author: Maciej Paprota
% Reference: M. Paprota. 2023. A twin wavemaker model for liquid sloshing in a rectangular tank. Ocean Engineering, 272, 113919

eta = eta_hat*cos(kappa'*x); % eq. (2)
end
