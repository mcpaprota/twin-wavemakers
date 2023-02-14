% Plotting linear sloshing in a 2D tank
% Author: Maciej Paprota
% Reference: M. Paprota. 2023. A twin wavemaker model for liquid sloshing in a rectangular tank. Ocean Engineering, 272, 113919

clear, clc
set(0,'defaulttextinterpreter','latex');
% Inialization:
g = 9.8145; % gravity acceleration (m/s^2)
l = 8; d = 2; % fluid domain size (m x m)
T = 3; H = 0.05; % wave period (s) and height (m)
nT = 20; nTi = 3; % number of wave periods and start-up periods (s)
Tdt = 1000; % number of time steps per wave period
K = 16; % number of eigenvalues/points in space
theta = pi/2; % wavemaker motion phase lag
[chi, xi, zeta, dt, nt] = paddleMotionRegular(d,H,T,nT,nTi,Tdt,g,theta);
[lambda, eta_hat, phi_hat] = twinWavemakerEulerMod(l,d,xi,zeta,dt,K,g);
% Plotting output:
figure('WindowState','maximized')
for n=1:Tdt/10:nt
    dx = l/40; x = chi(n):dx:l+chi(n); % fluid domain
    dz = d/20; z = (-d:dz:d)'; % fluid domain
    % free-surface elevation
    eta = freeSurfaceElevation(eta_hat(n,:),lambda,x);
    Eta = repmat(eta,length(z),1);
    % velocities
    [U, W, X, Z] = velocityComponents(l,d,xi(n),phi_hat(n,:),lambda,x,z,0);
    % removing velocities above free surface
    U(Z>Eta) = NaN; W(Z>Eta) = NaN;
    % water
    area(x,eta,-d,'FaceColor', [162 229 235]/255)
    axis([-abs(max(chi)) l+abs(max(chi)) -d d])
    hold on
    % velocity arrows (normalized by \sqrt(gd))
    quiver(X,Z,U/sqrt(g*d),W/sqrt(g*d),'AutoScale','off')
    % tank walls
    plot([chi(n) chi(n)],[-d d/2],'k','LineWidth',5)
    plot([l+chi(n) l+chi(n)],[-d d/2],'k','LineWidth',5)
    % tank bottom
    plot([chi(n) l+chi(n)],[-d -d],'k','LineWidth',5)
    set(gca,'Clipping','off')
    hold off
    daspect([1 1 1])
    title(['$t=$' num2str(dt*(n-1)) ' s'])
    xlabel('$x$ (m)'), ylabel('$z$ (m)')
    drawnow
end
