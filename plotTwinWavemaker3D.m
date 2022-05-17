% Plotting linear sloshing in a 3D tank
% Author: Maciej Paprota
% Reference: M. Paprota. 2022. A twin wavemaker model for liquid sloshing in a rectangular tank.

clear, clc
set(0,'defaulttextinterpreter','latex');
% Inialization:
g = 9.8145; % gravity acceleration (m/s^2)
l_1 = 8;l_2 = 8; d = 2; % fluid domain size (m x m)
T = 2; H = 0.2; % wave period (s) and height (m)
nT = 20; nTi = 3; % number of wave periods and start-up periods (s)
Tdt = 1000; % number of time steps per wave period
K = 16; % number of eigenvalues/points in space
theta_1 = 0; theta_2 = pi/2;
[chi_1, xi_1, zeta_1, dt, nt] = paddleMotionRegular(d,H,T,nT,nTi,Tdt,g,theta_1);
[lambda_1, eta_hat_1, phi_hat_1] = twinWavemakerEulerMod(l_1,d,xi_1,zeta_1,dt,K,g);
[chi_2, xi_2, zeta_2] = paddleMotionRegular(d,H,T,nT,nTi,Tdt,g,theta_2);
[lambda_2, eta_hat_2, phi_hat_2] = twinWavemakerEulerMod(l_2,d,xi_2,zeta_2,dt,K,g);
% Plotting output:
figure('WindowState','maximized')
colormap bone
d_w = 0.3;
for n=1:Tdt/10:nt
    dx = l_1/20; x = chi_1(n):dx:l_1+chi_1(n); % fluid domain
    dy = l_2/20; y = chi_2(n):dy:l_2+chi_2(n); % fluid domain
    dz = d/5; z = (-d:dz:d)'; % fluid domain
    % free-surface elevation
    eta_1 = freeSurfaceElevation(eta_hat_1(n,:),lambda_1,x);
    eta_2 = freeSurfaceElevation(eta_hat_2(n,:),lambda_2,y);
    eta = repmat(eta_1,length(y),1)+repmat(eta_2,length(x),1)';
    Eta = repmat(eta,1,1,length(z));
    % velocities
    [U, W_1, X, Z] = velocityComponents(l_1,d,xi_1(n),phi_hat_1(n,:),lambda_1,x,z,0);
    [V, W_2, Y] = velocityComponents(l_2,d,xi_2(n),phi_hat_2(n,:),lambda_2,y,z,0);
    % removing velocities above free surface
    U = permute(repmat(U,1,1,length(y)),[3 2 1]);
    V = permute(repmat(V,1,1,length(x)),[2 3 1]);
    W = permute(repmat(W_1,1,1,length(y)),[3 2 1])+...
        permute(repmat(W_2,1,1,length(x)),[2 3 1]);
    X = permute(repmat(X,1,1,length(y)),[3 2 1]);
    Y = permute(repmat(Y,1,1,length(x)),[2 3 1]);
    Z = permute(repmat(Z,1,1,length(y)),[3 2 1]);
    U(Z>Eta) = NaN; W(Z>Eta) = NaN;
    % water
    surf(X(:,:,1),Y(:,:,1),eta,'FaceLighting','gouraud',...
        'FaceColor',[100/255 200/255 225/255],'FaceAlpha',0.2,...
        'EdgeColor',[0.5 0.5 0.5])
    axis([-abs(max(chi_1)) l_1+abs(max(chi_1)) ...
        -abs(max(chi_2)) l_2+abs(max(chi_2)) -d d])
    hold on
%     % velocity arrows (normalized by \sqrt(gd))
%     quiver(X,Z,U/sqrt(g*d),W/sqrt(g*d),'AutoScale','off')
    % tank wall
    [vert, fac] = tankCoords([l_1/2+chi_1(n) l_2/2+chi_2(n) -d/4],l_1,l_2,1.5*d);
    patch('Vertices', vert, 'Faces', fac,...
        'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.3);
    quiver3(X,Y,Z,U/sqrt(g*d),V/sqrt(g*d),W/sqrt(g*d),'AutoScale','off',...
        'Color','red')
    set(gca,'Clipping','off')
    hold off
    daspect([1 1 1])
    title(['$t=$' num2str(dt*(n-1)) ' s'])
    xlabel('$x$ (m)'), ylabel('$y$ (m)'), zlabel('$z$ (m)')
    drawnow
end