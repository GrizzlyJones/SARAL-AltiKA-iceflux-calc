function plotLeads(Xq, Yq, pP_class)
warning('Implamentation of grid data will change in future versions');
figure;
title('Pulse Peakniss', 'FontSize', 18);
hold on
m_pcolor(Xq, Yq, pP_class);
shading flat;
caxis([0, 1]);
colorbar('Ticks', [0,1], 'TickLabels', {'No lead', 'Lead'})
colormap(parula(2));
m_grid;
end