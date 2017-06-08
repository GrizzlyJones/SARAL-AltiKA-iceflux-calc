function plotFlow(fluxgate, flow)

figure;
plot(fluxgate.profile.cumStep/1e3, flow/1e3);
title('Flow Distribution', 'FontSize', 18);
xlabel('Distance along fluxgate [km]');
ylabel('Flow rate [km day^{-1}]');

end