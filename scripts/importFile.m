addpath(fullfile(matlabroot, 'toolbox', 'matlab', 'm_map'));

cd('C:\Users\JonathanMemko\OneDrive\Dokumenter\DTU\6. semester\Bachelor\ALTIKA\2016');

lis = dir();
lis(1:2) = [];

LON = [-25, 20];
LAT = [76, 82];

m_proj('albers equal-area', 'long', LON, 'lat', LAT, 'rectbox', 'off');

lon = [];
lat = [];
ssha = [];

[Xq, Yq] = meshgrid(LON(1):.5:LON(2), LAT(1):.01:LAT(2));

for j = 1:1
    
    tmpLon = ncread(lis(j).name,'lon');
%     tmpLon = tmpLon * ncreadatt(lis(i).name, 'lon', 'scale_factor');
    tmpLon = rem((tmpLon + 180), 360) - 180;
    filLon = -25 < tmpLon & tmpLon < 20;
    tmpLat = ncread(lis(j).name,'lat');
%     tmpLat = tmpLat * ncreadatt(lis(i).name, 'lat', 'scale_factor');
    filLat = 76 < tmpLat & tmpLat < 82;
    tmpSsha = ncread(lis(j).name, 'waveforms_40hz');
%     tmpSsha = tmpSsha * ncreadatt(lis(i).name, 'tracker_40hz', 'scale_factor');
%     tmpSsha = tmpSsha + ncreadatt(lis(i).name, 'tracker_40hz', 'add_offset');

    filter = filLat & filLon;
    
    lon = [lon; tmpLon(filter)];
    lat = [lat; tmpLat(filter)];
    ssha = cat(3, ssha, tmpSsha(:,:,filter));
end

% Vq = griddata(lon, lat, ssha, Xq, Yq);
% 
% figure;
% hold on
% m_pcolor(Xq, Yq, Vq);
% % m_scatter(lon, lat, 5, ssha);
% shading flat;
% m_gshhs('lc', 'color', 'k');
% m_grid;

% figure;
% m_scatter(lon, lat, 5, ssha);
% shading flat;
% m_gshhs('lc', 'color', 'k');
% m_grid;

% figure('units', 'normalized', 'outerposition', [0 0 1 1]);
% % figure
% k = 1;
% for i = 1:5
%     for j = 1:length(ssha(1,:,1))
%         plot(ssha(:,j,i));
%         title(strcat(int2str(i), ', ', int2str(j)));
%         if sum(ssha(:,j,i)) > 0
%             [M, I] = max(ssha(:,j,i));
%             foo = find(ssha(1:I,j,i) < M/2);
%             N = foo(end);
%             n = N + (M/2 - ssha(N, j,i)) / (ssha(N+1, j, i) - ssha(N,j,i));
%             line([n n], get(gca, 'YLim'), 'color', 'g', 'LineStyle', '-');
%         end
% %         hax = axes;
%         xlim([1,128]);
%         line([64 64], get(gca, 'YLim'), 'color', 'r', 'LineStyle', '--');
% %         drawnow
%         F(k) = getframe;
%         k = k +1;
%     end
% end
% close all;
figure;
i = 23;
j = 53;
plot(ssha(:,i,j));
[M, I] = max(ssha(:,i,j));
foo = find(ssha(1:I,i,j) < M/2);
N = foo(end);
n = N + (M/2 - ssha(N, i, j)) / (ssha(N+1, i, j) - ssha(N,i,j));
line([n n], get(gca, 'YLim'), 'color', 'b', 'LineStyle', '-');
line([64 64], get(gca, 'YLim'), 'color', 'r', 'LineStyle', '--');