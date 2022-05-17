function [u, w, X, Z] = ...
    velocityComponents(l,d,xi,phi_hat,kappa,x,z,notTwin)
% calculates component velocities as spatial derivatives of potential function based on velocity of a paddle(s) and amplitudes and eigenvalues of cosine expansion.
%
% Input data:
% xi - paddle(s) velocity (m/s)
% phi_hat - solution amplitudes (m)
% kappa(lambda) - solution aigenvalues (1/m)
% x - longitudinal coordinates (m)
% z - vertical coordinates (m)
% notTwin - negative twin wavemaker indicator 0 - false, 1 - true
%
% Output data:
% u, w - horizontal and vertical velocity components (m/s)
%
% Author: Maciej Paprota
% Reference: M. Paprota. 2022. A twin wavemaker model for liquid sloshing in a rectangular tank.
%
% Initialization
[X,Z,Kappa] = meshgrid(x,z,kappa);
[~,~,Phi_hat] = meshgrid(x,z,phi_hat);
u = -notTwin*xi*X(:,:,1)/l+xi-... % conditionally resolved eqs. (2.20) and (2.35)
    sum(Phi_hat.*Kappa.*cosh(Kappa.*(Z+d))./...
    cosh(Kappa*d).*sin(Kappa.*X),3);
w = notTwin*xi*(Z(:,:,1)+d)/l+sum(Phi_hat.*Kappa.*sinh(Kappa.*(Z+d))./... % conditionally resolved eqs. (2.21) and (2.36)
    cosh(Kappa*d).*cos(Kappa.*X),3);
X = X(:,:,1); Z = Z(:,:,1);
end