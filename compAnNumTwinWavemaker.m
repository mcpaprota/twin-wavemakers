% Comparison of analytical and numerical solution to linear twin wavemaker problem
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
T = 2*pi/sigma; % wavemaker period (s)
dt = T/Tdt; % time increment (s)
nT = 10; % number of periods for analysis
t = (0:dt:T*nT)'; % time domain vector (s)
% wavemaker paddle kinematics
chi = chi_hat*sin(sigma*t); % paddle displacement (m)
xi = chi_hat*sigma*cos(sigma*t); % paddle velocity (m/s)
zeta = -chi_hat*sigma^2*sin(sigma*t); % paddle acceleration (m^2/s)
% analytical solution
[lambdaA, eta_hatA, phi_hatA] = ...
    twinWavemakerAnalytical(l,d,chi_hat,sigma,t,K,g);
% numerical solution
[lambdaN, eta_hatN, phi_hatN] = twinWavemakerEulerMod(l,d,xi,zeta,dt,K,g);
% plotting solution amplitudes
subplot(2,1,1)
plot(lambdaA,eta_hatA(end,:),'k',lambdaN,eta_hatN(end,:),'b')
xlabel('$\kappa$'), ylabel('$\hat{\eta}$ (m)')
subplot(2,1,2)
plot(lambdaA,phi_hatA(end,:),'k',lambdaN,phi_hatN(end,:),'b')
xlabel('$\kappa$'), ylabel('$\hat{\varphi}$ (m$^2$/s)')
