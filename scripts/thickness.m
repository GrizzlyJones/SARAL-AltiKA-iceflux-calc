function [iceSheet] = thickness( iceSheet )
%THICKNESS Calculates icethickness from freeboard
%   h_i = thickness(H_SF) calculates ice thickness assuming altimitry is
%   done using laser. The equation is based on hydrostatic balance. H_SF is
%   the freeboard height.
% 
%   Source: "Thickness and density of snow-covered sea ice and hydrostatic
%   equilibrium assumption from in situ measurements in Fram Straint, the
%   Barents Sea and the Svalbard coast"
%       by Forsström, Gerland and Pedersen

 
rho_w = 1020;   % Density of water in kg/m^3 
rho_i = 910;    % Density of ice in kg/m^3
rho_s = 300;    % Density of snow in kg/m^3
h_s = 0.147;    % Snow height from ESA CryoSat 2


for i = 1:length(iceSheet)
    h_sf = iceSheet(i).freeboard;
    iceSheet(i).iceThickness = h_sf .* rho_w/(rho_w - rho_i) + h_s * (rho_s-rho_w)/(rho_w-rho_i);
end

end