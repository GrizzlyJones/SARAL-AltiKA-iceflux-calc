function freeboard_line = cryoSat(filePath, var, fluxgate, LON, LAT)

% Tmp import of Lon and Lat
tmpLon = ncread(filePath, 'longitude');
tmpLon = rem((tmpLon + 180), 360) - 180;
filLon = (LON(1) - 1) < tmpLon & tmpLon < (LON(2) + 1);

tmpLat = ncread(filePath,'latitude');
filLat = (LAT(1) - 1) < tmpLat & tmpLat < (LAT(2) + 1);

% Filter for earth
filter = filLat & filLon;

tmpFreeboard = ncread(filePath, var);

lon = tmpLon(filter);
lat = tmpLat(filter);
freeboard = tmpFreeboard(filter);

%%
freeboard_q = griddata(double(lon), double(lat), double(freeboard), Xq, Yq);

freeboard_line = interp2(Xq, Yq, freeboard_q, fluxgate.profile.lon, fluxgate.profile.lat);

% %%
% figure;
% subplot(2,1,1);
% title('Retracked, scatter');
% m_scatter(lon, lat, 10, freeboard, 'filled');
% colorbar;
% m_grid;
% subplot(2,1,2);
% title('Retracked SLA');
% m_pcolor(Xq, Yq, freeboard_q);
% shading flat;
% colorbar;
% m_grid;
% 
% %%
% figure;
% plot(fluxgate.profile.cumStep/1e3, freeboard_line);
end