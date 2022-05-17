function kappa = dispersionRelation(d,omega,g)
% calculates wave number based on linear dispertion relationship
%
% Input data:
% d - water depth (m)
% omega - vector of angular frequencies (1/s)
% g - acceleration due to gravity (m/s^2)
%
% Output data:
% kappa - wave numbers (1/m)
%
% Author: Maciej Paprota
% Reference: M. Paprota. 2022. A twin wavemaker model for liquid sloshing in a rectangular tank

eps = 10^(-12); % tolerance
kappa = omega.*omega/g; % initial kappa values
% Iteration loop:
while max(abs(kappa - omega.*omega/g./tanh(kappa*d)))>eps
    kappa = omega.*omega/g./tanh(kappa*d);
end
end