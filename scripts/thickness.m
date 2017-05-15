function [T] = thickness(freeboard)

% Sea ice thickness


T = freeboard .* rho_w/(rho_w - rho_i) + Z * (rho_s-rho_w)/(rho_w-rho_i) + P * rho_w/(rho_w-rho_i);

% IceBridge L4 Sea Ice Freeboard, Snow Depth, and Thickness giver nogle
% andre værdier for nedenstående konstanter. 
rho_w = 1020;   % Density of water in kg/m^3 
rho_i = 910;    % Density of ice in kg/m^3
rho_s = 300;    % Density of snow in kg/m^3
Z = 0.3;        % Snow height from ESA cryosat2
P = 0.15;       % Penetration from ESA cryosat2


end