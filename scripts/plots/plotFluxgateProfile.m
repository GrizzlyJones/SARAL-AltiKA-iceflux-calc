function plotFluxgateProfile(Xq, Yq, ssha_q, fluxgate, varargin)
%PLOTFLUXGATEPROFILE plots the fluxgate and freeboards
%   void = plotFluxgate(XQ, YQ, SSHA_Q, FLUXGATE, VARARGIN) plots the
%   FLUXGATE profile on a map alongside SSHA_Q
% 
%   See also INITFLUXGATE, INTERPPROFILE, FREEBOARDANALYSIS

warning('Implamentation of grid data will change in future versions');

profile = fluxgate.profile;

figure;
m_pcolor(Xq, Yq, ssha_q);
shading flat;
cb = colorbar;
cb.Label.String = 'Height [m]';
m_grid;
m_track(profile.lon, profile.lat);
title('Retracked SLA');

% if ~isempty(varargin)
%     fnam = sprintf('figures/fluxgate_profile_%s', varargin{2});
%     hgexport(gcf, fnam, varargin{1});
% end
end

