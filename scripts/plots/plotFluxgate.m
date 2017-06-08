function plotFluxgate(fluxgate,  freeboard)
%PLOTFLUXGATE plots the fluxgate and freeboards
%   void = plotFluxgate(XQ, YQ, SSHA_Q, FLUXGATE, FREEBOARD) plots the
%   FLUXGATE on a map alongside SSHA_Q, as well as plotting the FREEBOARD
% 
%   See also INITFLUXGATE, INTERPPROFILE, FREEBOARDANALYSIS

data = fluxgate.data;

figure;
dist = fluxgate.profile.cumStep/10e3;

subplot(2,1,1);
hold on
class = data.class;
class(class == 0) = nan;
class(~isnan(class)) = -1;
plot(dist, class, 'color', 'r', 'marker', '*', 'lineWidth', 1, 'lineStyle', 'none');
plot(dist, data.sla, 'color', [0 0.4470 0.7410]);
xlim([0, dist(end)]);
xlabel('Distance [km]');
ylabel('Height [m]');
title('Fluxgate SLA', 'fontSize', 18);
legend('Leads', 'SLA', 'location', 'southeast');

subplot(2,1,2);
hold on
for i = 1:length(freeboard)
    plot(dist(freeboard(i).index.start:freeboard(i).index.stop), ...
         freeboard(i).freeboard);
end
xlim([0, dist(end)]);
xlabel('Distance [km]');
ylabel('Height [m]');
title('Fluxgate Freeboards', 'fontSize', 18);
end

