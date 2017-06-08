function plotSLA(Xq, Yq, ssha_q, title_)
warning('Implamentation of grid data will change in future versions');

figure;
m_pcolor(Xq, Yq, ssha_q);
cb = colorbar;
cb.Label.String = 'Height [m]';
shading flat;
m_grid;
title(title_, 'fontSize', 18);
end