function plotArea(LON, LAT, fluxgate)
%PLOTFLUXGATEPROFILE plots the fluxgate and freeboards
%   void = plotFluxgate(XQ, YQ, SSHA_Q, FLUXGATE, VARARGIN) plots the
%   FLUXGATE profile on a map alongside SSHA_Q
% 
%   See also INITFLUXGATE, INTERPPROFILE, FREEBOARDANALYSIS

m_proj('stereographic', 'long', 0, 'lat', 80, 'radius', 10);
% m_proj('satellite', 'long', 0, 'lat', 80, 'alt', .1);
% m_proj('orthographic', 'long', 0, 'lat', 80);
% m_proj('Azimuthal Equal-Area', 'long', 0, 'lat', 80, 'rad', 10);
% m_proj('Mercator', 'long', [-15 15], 'lat', [75 85]);
% m_proj('Mollweide', 'long', [-15 15], 'lat', [75 85]);

figure;
m_gshhs('lc', 'patch', [0.9 0.9 1], 'edgeColor', [0.3 0.3 0.3]);
title('Fram Strait', 'FontSize', 18);
bndry_lon = [LON(1), LON(2), LON(2) LON(1), LON(1)];
bndry_lat = [LAT(1), LAT(1), LAT(2), LAT(2), LAT(1)];
m_line(bndry_lon, bndry_lat, 'color', 'k', 'linewi', 1);
m_hatch(bndry_lon, bndry_lat, 'single', 45, 5, 'color', [.6,.6,.6]);

m_track(fluxgate.profile.lon, fluxgate.profile.lat, 'color', 'r');
hold on
foo = plot(nan, nan, 'color', 'r');
legend(foo, 'Fluxgate')
m_grid('ytick', 5);

end

