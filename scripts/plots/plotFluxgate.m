function plotFluxgate(Xq, Yq, ssha_q, fluxgate,  freeboard)
%PLOTFLUXGATE Summary of this function goes here
%   Detailed explanation goes here
figure;
subplot(1,2,1);
m_pcolor(Xq, Yq, ssha_q);
shading flat;
colorbar;
m_grid;
m_track(fluxgate.lon, fluxgate.lat);
title('Retracked SLA');

subplot(2,2,2);
hold on
plot(fluxgate.sla);
plot(fluxgate.class/4, 'r');
title('From cubic interpolation');

subplot(2,2,4);
hold on
for i = 1:length(freeboard)
    plot(freeboard(i).index.start:freeboard(i).index.stop, ...
         freeboard(i).freeboard);
end
xlim([0 fluxgate.steps]);
title('Freeboards');
end

