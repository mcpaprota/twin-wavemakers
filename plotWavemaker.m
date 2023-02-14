% Plotting mechanically generated waves in a 2D wave flume
% Author: Maciej Paprota
% Reference: M. Paprota. 2023. A twin wavemaker model for liquid sloshing in a rectangular tank. Ocean Engineering, 272, 113919

clear, clc
set(0,'defaulttextinterpreter','latex');
% Inialization:
g = 9.8145; % gravity acceleration (m/s^2)
l = 80; d = 2; % fluid domain size (m x m)
T = 2; H = 0.2; % wave period (s) and height (m)
nT = 20; nTi = 3; % number of wave periods and start-up periods (s)
Tdt = 1000; % number of time steps per wave period
I = 128; % number of eigenvalues/points in space
theta = 0; % wavemaker motion phase lag
[chi, xi, zeta, dt, nt] = paddleMotionRegular(d,H,T,nT,nTi,Tdt,g,theta);
[kappa, eta_hat, phi_hat] = wavemakerEulerMod(l,d,xi,zeta,dt,I,g);
% Plotting output:
figure('WindowState','maximized')
for n=1:Tdt/10:nt
    dx = (l-chi(n))/80; x = chi(n):dx:l; % fluid domain
    dz = d/20; z = (-d:dz:d)'; % fluid domain
    % free-surface elevation
    eta = freeSurfaceElevation(eta_hat(n,:),[0 kappa],x);
    Eta = repmat(eta,length(z),1);
    % velocities
    [U, W, X, Z] = velocityComponents(l,d,xi(n),phi_hat(n,:),kappa,x,z,1);
    % removing velocities above free surface
    U(Z>Eta) = NaN; W(Z>Eta) = NaN;
    % water
    area(x,eta,-d,'FaceColor', [162 229 235]/255)
    axis([-abs(max(chi)) l -d d])
    hold on
    % velocity arrows (normalized by \sqrt(gd))
    quiver(X,Z,U/sqrt(g*d),W/sqrt(g*d),'AutoScale','off')
    % wavemaker paddle
    plot([chi(n) chi(n)],[-d d/2],'k','LineWidth',5)
    % far end wall
    plot([l l],[-d d/2],'k','LineWidth',5)
    % bottom
    plot([-abs(max(chi)) 1.002*l],[-d -d],'k','LineWidth',5)
    set(gca,'Clipping','off')
    hold off
    daspect([1 1 1])
    title(['$t=$' num2str(dt*(n-1)) ' s'])
    xlabel('$x$ (m)'), ylabel('$z$ (m)')
    drawnow
end
