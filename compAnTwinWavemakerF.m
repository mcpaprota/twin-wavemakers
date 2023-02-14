% Comparison of twin wavemaker problem solution and Faltinsen (1978) solution to sloshing
% Author: Maciej Paprota
% Reference: M. Paprota. 2023. A twin wavemaker model for liquid sloshing in a rectangular tank. Ocean Engineering, 272, 113919

clc, clear
set(0,'defaulttextinterpreter','latex');
% Initialization:
g = 9.8145; % gravity acceleration (m/s^2)
l = 10; d = 2; % fluid domain size (m x m)
chi_hat = 0.05; % wavemaker displacement amplitude (m)
R = 2.5; % mode number (integer values generate resonance)
sigma = sqrt(g*(2*R-1)*pi/l*tanh((2*R-1)*pi/l*d)); % wavemaker frequency (rad/s)
Tdt = 1000; % number of time steps per wave period
K = 20; % number of eigenvalues
dx = l/200; x = 0:dx:l; % fluid domain (m)
x_prime = x-l/2; % Faltinsen (1978) fluid domain
T = 2*pi/sigma; % wavemaker period (s)
dt = T/Tdt; % time increment (s)
nT = 1; % number of periods for analysis
t = (0:dt:T*nT)'; % time domain vector (s)
% wavemaker paddle kinematics
chi = chi_hat*sin(sigma*t); % paddle displacement (m)
xi = chi_hat*sigma*cos(sigma*t); % paddle velocity (m/s)
zeta = -chi_hat*sigma^2*sin(sigma*t); % paddle acceleration (m^2/s)
% analytical solution
[~, ~, phi_hatA] = ...
    twinWavemakerAnalytical(l,d,chi_hat,sigma,t,K,g);
[lambdaA, phi_hatAF] = ...
    twinWavemakerAnalyticalF(l,d,chi_hat,sigma,t,K,g);
% plotting velocity potential function
for n=1:length(t)
    phiA = chi_hat*sigma*cos(sigma*t(n))*(x-l/2)+phi_hatA(n,:)*cos(lambdaA'*x);
    phiAF = chi_hat*sigma*cos(sigma*t(n))*(x_prime)+phi_hatAF(n,:)*sin(lambdaA'*(x_prime));
    plot(x,phiA(end,:),'k',x,phiAF(end,:),'b')
    xlabel('$x$'), ylabel('$\varphi$ (m$^2$/s)')
    drawnow
end
