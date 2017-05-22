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

dist = fluxgate.profile.cumStep/10e3;

subplot(2,2,2);
hold on
plot(dist, data.sla);
plot(dist, data.class/4, 'r');
xlim([0, dist(end)]);
xlabel('Distance [km]');
ylabel('Height [m]');
title('From cubic interpolation');

subplot(2,2,4);
hold on
for i = 1:length(freeboard)
    plot(dist(freeboard(i).index.start:freeboard(i).index.stop), ...
         freeboard(i).freeboard);
end
xlim([0, dist(end)]);
xlabel('Distance [km]');
ylabel('Height [m]');
title('Freeboards');
end

