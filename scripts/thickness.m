function [iceSheet] = thickness( iceSheet , type )
%THICKNESS Calculates icethickness from freeboard
%   h_i = thickness(ICESHEET, TYPE) calculates ice thickness assuming altimitry is
%   done using laser. The equation is based on hydrostatic balance. h_sf is
%   the freeboard height found within iceSheet. Either laser or radar
%   algorithme can be used, defaults to laser.
% 
%   Source: "Thickness and density of snow-covered sea ice and hydrostatic
%   equilibrium assumption from in situ measurements in Fram Straint, the
%   Barents Sea and the Svalbard coast"
%       by Forsström, Gerland and Pedersen

if nargin < 2
    type = 'laser';
end
 
rho_w = 1020;   % Density of water in kg/m^3 
rho_i = 910;    % Density of ice in kg/m^3
rho_s = 300;    % Density of snow in kg/m^3
h_s = 0.147;    % Snow height from ESA CryoSat 2

switch lower(type)
    case 'laser'
        for i = 1:length(iceSheet)
            h_sf = iceSheet(i).freeboard;
            iceSheet(i).iceThickness = h_sf .* rho_w/(rho_w - rho_i) + h_s * (rho_s-rho_w)/(rho_w-rho_i);
        end
    case 'radar'
        for i = 1:length(iceSheet)
            h_sf = iceSheet(i).freeboard;
            iceSheet(i).iceThickness = h_sf .* rho_w/(rho_w - rho_i) + h_s * (rho_s)/(rho_w-rho_i);
        end
    otherwise
        error('Error. \nType must be ''radar'' or ''laser'', note %s', type);
end
end