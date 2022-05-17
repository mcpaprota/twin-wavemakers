function [vert,fac] = tankCoords(C, X, Y, Z)
% defines coordinates of a patch object according to vertices and faces of a tank
%
% Input data:
% C - center position (m)
% X, Y, Z - paddle lengths along x, y, z axes (m)
% 
% Output data:
% vert - vertices coordinates (m)
% fac - faces indicators
%
% Author: Maciej Paprota
% Reference: M. Paprota. 2022. A twin wavemaker model for liquid sloshing in a rectangular tank

x_vert = [0 1 1 0 0 1 1 0]*X-X/2+C(1);
y_vert = [0 0 1 1 0 0 1 1]*Y-Y/2+C(2);
z_vert = [0 0 0 0 1 1 1 1]*Z-Z/2+C(3);
vert = [x_vert' y_vert' z_vert'];
fac = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4];
end