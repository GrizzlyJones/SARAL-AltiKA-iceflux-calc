%% Init
clc; close all; clear;

%% File mangagement
addpath(fullfile(matlabroot, 'toolbox', 'matlab', 'm_map')); % m_maps
addpath(fullfile(pwd,'scripts'));                            % used scripts
altikaFiles = fullfile(pwd,'ALTIKA','2016');                 % data

% All data files
lis = dir(altikaFiles);
lis(1:2) = [];

%% Fram Strait
% LON = [-25, 20];
LON = [-10, 10];
LAT = [76, 82];

m_proj('albers equal-area', 'long', LON, 'lat', LAT, 'rectbox', 'off');

lon = [[], []];
lat = [[], []];
wave = [[], []];
agc = [[], []];
track = [[], []];
alt = [[], []];

modeled_instr_corr = [[], []];
doppler_corr = [[], []];

model_dry_tropo_corr = [[], []];
model_wet_tropo_corr = [[], []];
iono_corr_gim = [[], []];
sea_state_bias = [[], []];

range = [[], []];

% [Xq, Yq] = meshgrid(LON(1):.5:LON(2), LAT(1):.01:LAT(2));

%% Load data
% for j = 1:10;
%     filePath = fullfile(altikaFiles, lis(j).name);
%     tmpLon = ncImport(filePath, 'lon_40hz');
%     tmpLon = rem((tmpLon + 180), 360) - 180;
%     filLon = LON(1) < tmpLon & tmpLon < LON(2);
%     tmpLat = ncImport(filePath,'lat_40hz');
%     filLat = LAT(1) < tmpLat & tmpLat < LAT(2);
%     
%     tmpWave = ncImport(filePath, 'waveforms_40hz');
%     tmpAGC = ncImport(filePath, 'agc_40hz');
%     tmpRange = ncImport(filePath, 'range_40hz');
%     tmpAlt = ncImport(filePath, 'alt_40hz');
%     
%     filAlt = tmpAlt < mean(tmpAlt(:));
%     filRange = tmpRange < mean(tmpRange(:));
%     
%     filter = filLat & filLon & filAlt & filRange;
%     
%     lon = [lon; tmpLon(filter)];
%     lat = [lat; tmpLat(filter)];
%     wave = [wave, tmpWave(:,filter)];
%     agc = [agc; tmpAGC(filter)];
%     range = [range; tmpRange(filter)];
%     alt = [alt; tmpAlt(filter)];
% end
for j = 3:4;
    filePath = fullfile(altikaFiles, lis(j).name);
    tmpLon = ncread(filePath, 'lon_40hz');
    tmpLon = rem((tmpLon + 180), 360) - 180;
    filLon = LON(1) < tmpLon & tmpLon < LON(2);
    tmpLat = ncread(filePath,'lat_40hz');
    filLat = LAT(1) < tmpLat & tmpLat < LAT(2);
    
    tmpWave = ncread(filePath, 'waveforms_40hz');
    tmpAGC = ncread(filePath, 'agc_40hz');
    tmpTrack = ncread(filePath, 'tracker_40hz');
    tmpAlt = ncread(filePath, 'alt_40hz');
    
    tmpModeled_intr_corr = ncread(filePath, 'modeled_instr_corr_range');
    tmpDoppler_corr = ncread(filePath, 'doppler_corr');
    
    tmpModel_dry_tropo_corr = ncread(filePath, 'model_dry_tropo_corr');
    tmpModel_wet_tropo_corr = ncread(filePath, 'model_wet_tropo_corr');
    tmpIono_corr_gim = ncread(filePath, 'iono_corr_gim');
    tmpSea_state_bias = ncread(filePath, 'sea_state_bias');
    
    tmpRange = ncread(filePath, 'range_40hz');
    
    filter = filLat & filLon;
    
    lon = [lon; tmpLon(filter)];
    lat = [lat; tmpLat(filter)];
    wave = [wave, tmpWave(:,filter)];
    agc = [agc; tmpAGC(filter)];
    track = [track; tmpTrack(filter)];
    alt = [alt; tmpAlt(filter)];
    
    modeled_instr_corr = [modeled_instr_corr; tmpModeled_intr_corr(filter(1, :))];
    doppler_corr = [doppler_corr; tmpDoppler_corr(filter(1,:))];
    
    model_dry_tropo_corr = [model_dry_tropo_corr; tmpModel_dry_tropo_corr(filter(1, :))];
    model_wet_tropo_corr = [model_wet_tropo_corr; tmpModel_wet_tropo_corr(filter(1, :))];
    iono_corr_gim = [iono_corr_gim; tmpIono_corr_gim(filter(1, :))];
    sea_state_bias = [sea_state_bias; tmpSea_state_bias(filter(1, :))];
    
    range = [range; tmpRange(filter)];
    disp(j)
end

%% Calculate height
C_ntp = 64;
B_spc = 0.31;

N = length(wave(1,:));
C_rtrk_ocog = zeros(N, 1);
C_rtrk_pp_cog = zeros(N, 1);
mp = zeros(N,1);
pP = zeros(N,1);

for i = 1:N
    C_rtrk_ocog(i) = waveformAnalysis(wave(:,i), 'OCOG');
    C_rtrk_pp_cog(i) = waveformAnalysis(wave(:,i), 'PP_COG');
    mp(i) = maxPower(wave(:,i), agc(i));
    pP(i) = pulsePeakness(wave(:,i));
end

type = zeros(N, 1);
type(pP > 30) = 4;

type2 = zeros(N, 1);
% type2(40 <= pP < 50) = 2;
% type2(50 <= pP < 70) = 3;
type2(70 <= pP) = 4;

x = linspace(1,N, length(modeled_instr_corr));

track_corr = modeled_instr_corr + doppler_corr;
track_corr_40hz = interp1(x, track_corr, 1:N)';

alt_corr = model_dry_tropo_corr + model_wet_tropo_corr + iono_corr_gim + sea_state_bias;
alt_corr_40hz = interp1(x, alt_corr, 1:N)';

epoch_ocog = (C_ntp - C_rtrk_ocog) * B_spc;
epoch_pp_cog = (C_ntp - C_rtrk_pp_cog) * B_spc;

ssh_ocog = alt + alt_corr_40hz - track + track_corr_40hz + epoch_ocog;
ssh_pp_cog = alt + alt_corr_40hz - track + track_corr_40hz + epoch_pp_cog;

% ssh_ocog = alt - track + epoch_ocog;
% ssh_pp_cog = alt - track + epoch_pp_cog;

%% Plot
% [Xq, Yq] = meshgrid(LON(1):0.1:LON(2), LAT(1):0.01:LAT(2));
% Vq = griddata(lon, lat, ssh_ocog, Xq, Yq);
% 
% figure;
% hold on
% m_pcolor(Xq, Yq, Vq);
% shading flat;
% % m_scatter(lon, lat, 5, ssha);
% m_gshhs('lc', 'color', 'k');
% m_grid;

figure;
% m_pcolor(Xq, Yq, Vq);
% c = flag(4);
% colormap(c);
m_scatter(lon, lat, 10, 'filled');
% shading flat;
m_gshhs('ic', 'color', 'k');
m_grid;

% figure;
% plot(mp);

% foo = zeros(N,1);
% for i = 1:N
%     foo(i) = pulsePeakness(wave(:,i));
% end

%%
% figure
% hold on;
% plot(ssh_ocog);
% plot(ssh_pp_cog);
% legend('OCOG', 'PP_{COG}');
% 
figure;
subplot(2,1,1);
hold on
xDis = (1:N) * 0.180;
plot(xDis, ssh_pp_cog);
% plot(xDis, mp);
% plot(xDis, pP);
% plot(get(gca, 'XLim'), [0,0], '--k');
% plot(get(gca, 'XLim'), [70, 70], '--k');
set(gca, 'xDir', 'reverse', 'YAxisLocation', 'right');

xlabel('Distance traversed [km]');
ylabel('Height / Power [m / dB]');

legend('Height', 'Maximum Power', 'Location', 'southeast');

subplot(2,1,2);
colormap(flipud(flag(3)));
m_scatter(lon, lat, 10, type, 'filled');
m_gshhs('lc', 'color', 'k');
m_grid;

%%
% figure;
% subplot(2,3,1);
% plot(wave(:,epoch_ocog <= -6));
% title('dis <= -6');
% subplot(2,3,2);
% plot(wave(:,-6 < epoch_ocog <= -3));
% title('-6 < dis <= -3');
% subplot(2,3,3);
% plot(wave(:, -3 < epoch_ocog <= 0));
% title('-3 < dis <= 0');
% subplot(2,3,4);
% plot(wave(:,0 < epoch_ocog <= 3));
% title('0 < dis <= 3');
% subplot(2,3,5);
% plot(wave(:,3 < epoch_ocog <= 6));
% title('3 < dis <= 6');
% subplot(2,3,6);
% plot(wave(:,epoch_ocog > 6));
% title('dis > 6');

%% Test af Primary Peak
% for i = 1:1;
%     i = i*2500;
%     figure;
%     [pp_start, pp_stop, Th_start, Th_stop, d_1] = primaryPeak(wave(:,i));
%     hold on;
%     plot(1:128, wave(:,i), 'r')
%     plot(pp_start:pp_stop, wave(pp_start:pp_stop,i), 'k')
%     title(i);
% end