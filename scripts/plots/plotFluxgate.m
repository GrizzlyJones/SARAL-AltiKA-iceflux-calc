function plotFluxgate(Xq, Yq, ssha_q, fluxgate,  freeboard)
%PLOTFLUXGATE plots the fluxgate and freeboards
%   void = plotFluxgate(XQ, YQ, SSHA_Q, FLUXGATE, FREEBOARD) plots the
%   FLUXGATE on a map alongside SSHA_Q, as well as plotting the FREEBOARD
% 
%   See also INITFLUXGATE, INTERPPROFILE, FREEBOARDANALYSIS

warning('Implamentation of grid data will change in future versions');

profile = fluxgate.profile;
data = fluxgate.data;

figure;
subplot(1,2,1);
m_pcolor(Xq, Yq, ssha_q);
shading flat;
colorbar;
m_grid;
m_track(profile.lon, profile.lat);
title('Retracked SLA');

subplot(2,2,2);
hold on
plot(data.sla);
plot(data.class/4, 'r');
title('From cubic interpolation');

subplot(2,2,4);
hold on
for i = 1:length(freeboard)
    plot(freeboard(i).index.start:freeboard(i).index.stop, ...
         freeboard(i).iceThickness);
end
xlim([0 profile.steps]);
title('Freeboards');
end

