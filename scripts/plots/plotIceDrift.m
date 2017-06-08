function plotIceDrift(Xq, Yq, sla_pp_cog_q, velocities)
warning('Implamentation of grid data will change in future versions');

figure;
hold on
m_pcolor(Xq, Yq, sla_pp_cog_q);
scale = 0;
m_quiver(velocities.lon, velocities.lat, velocities.x/100000, velocities.y/100000, scale);
m_quiver(4.5, 79.25, .5, 0, scale,'color', 'k');
m_text(6, 79.13, '50 km day^{-1}', 'FontSize', 12, 'horizontalAlignment', 'center');
shading flat;
cb = colorbar;
cb.Label.String = 'Height [m]';
m_grid;
title('Ice drift', 'fontSize', 18);
end