function plotCryoSat(fluxgate, altikaH, cryosatH, title_)
figure;
hold on
% subplot(3,1,1)
plot(fluxgate.profile.cumStep/1e3, cryosatH, 'lineWidth', 1, 'color', 'r');
% title(sprintf('CryoSat 2 %s', strVar), 'fontSize', 18);
% xlabel('Distance along fluxgate [km]');
% ylabel('Height [m]');

% subplot(3,1,2)
plot(fluxgate.profile.cumStep/1e3, altikaH, 'lineWidth', 1, 'color', 'b');
% title(sprintf('SARAL/AltiKa %s', strVar), 'fontSize', 18);
% xlabel('Distance along fluxgate [km]');
% ylabel('Height [m]');

% subplot(3,1,3)
plot(fluxgate.profile.cumStep/1e3, altikaH - cryosatH, 'lineWidth', 1, 'color', [1 0 1]);
title(title_, 'fontSize', 18);
xlabel('Distance along fluxgate [km]');
ylabel('Height [m]');
legend('CryoSat', 'SARAL/AltiKa', 'Diference');
end