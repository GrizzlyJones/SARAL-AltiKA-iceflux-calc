%% Init
clc; close all; clear;

%% File mangagement
addpath(fullfile(matlabroot, 'toolbox', 'matlab', 'm_map')); % m_maps
addpath(fullfile(pwd,'scripts'));                            % used scripts
altikaFiles = 'E:\Altika';                                   % data
load('scripts/style');                                       % export type

%% Fram Strait
LON = [-10, 10];
% LAT = [76, 82];
LAT = [79, 82];

% Settings for map projection
m_proj('albers equal-area', 'long', LON, 'lat', LAT, 'rectbox', 'off');

% Init. data arrays
lon = [];
lat = [];
wave = [];
agc = [];
tracker = [];
alt = [];

modeled_instr_corr = [];
doppler_corr = [];

model_dry_tropo_corr = [];
rad_wet_tropo_corr = [];
iono_corr_gim = [];
sea_state_bias = [];

range = [];

mss = [];
ssha = [];

solidEarthTideHeight = [];
oceanTide = [];
poleTide = [];
invBarCorr = [];
HF = [];

n = 0;
N = [];

%% Load data
cycle = 32;
j = 109;

cycleName = sprintf('cycle_%03d', cycle);
cycleFile = fullfile(pwd,'data', strcat(cycleName, '.mat'));
cycleFilePath = fullfile(altikaFiles, cycleName);
lis = dir(cycleFilePath);
lis(1:2) = [];
% Extraction of the filepath
filePath = fullfile(cycleFilePath, lis(j).name);

% Tmp import of Lon and Lat
tmpLon = ncread(filePath, 'lon_40hz');
tmpLon = rem((tmpLon + 180), 360) - 180;
filLon = (LON(1) - 1) < tmpLon & tmpLon < (LON(2) + 1);
tmpLat = ncread(filePath,'lat_40hz');
filLat = (LAT(1) - 1) < tmpLat & tmpLat < (LAT(2) + 1);

% Filter for earth
filter = filLat & filLon;

% Skips iteration if no useful data is detected
if ~(any(filter(:))) || sum(filter(:)) < 5
    fprintf('%d skipped, no data\n', j);
    return;
end

% Tmp import of data variables
tmpWave = ncread(filePath, 'waveforms_40hz');
tmpAGC = ncread(filePath, 'agc_40hz');
tmpTracker = ncread(filePath, 'tracker_40hz');
tmpAlt = ncread(filePath, 'alt_40hz');

tmpModeled_intr_corr = ncread(filePath, 'modeled_instr_corr_range');
tmpDoppler_corr = ncread(filePath, 'doppler_corr');

tmpModel_dry_tropo_corr = ncread(filePath, 'model_dry_tropo_corr');
tmpRad_wet_tropo_corr = ncread(filePath, 'rad_wet_tropo_corr');
tmpIono_corr_gim = ncread(filePath, 'iono_corr_gim');
tmpSea_state_bias = ncread(filePath, 'sea_state_bias');

tmpRange = ncread(filePath, 'range_40hz');

tmpMss = ncread(filePath, 'mean_sea_surface');
tmpSsha = ncread(filePath, 'ssha');

tmpSolidEarthTideHeight = ncread(filePath, 'solid_earth_tide');
tmpOceanTide = ncread(filePath, 'ocean_tide_sol2');
tmpPoleTide = ncread(filePath, 'pole_tide');
tmpInvBarCorr = ncread(filePath, 'inv_bar_corr');
tmpHF = ncread(filePath, 'hf_fluctuations_corr');

tmpN = length(tmpLon(filter));
x = linspace(1,length(tmpHF(filter(1,:))), tmpN);

% Filtration of found data, saving for later use
lon = cat(3, lon, tmpLon(filter));
lat = cat(3, lat, tmpLat(filter));
wave = cat(3, wave, tmpWave(:,filter));
agc = cat(3, agc, tmpAGC(filter));
tracker = cat(3, tracker, tmpTracker(filter));
alt = cat(3, alt, tmpAlt(filter));

modeled_instr_corr = cat(3, modeled_instr_corr, interp1(tmpModeled_intr_corr(filter(1, :)), x)');
doppler_corr = cat(3, doppler_corr, interp1(tmpDoppler_corr(filter(1,:)), x)');

model_dry_tropo_corr = cat(3, model_dry_tropo_corr, interp1(tmpModel_dry_tropo_corr(filter(1, :)), x)');
rad_wet_tropo_corr = cat(3, rad_wet_tropo_corr, interp1(tmpRad_wet_tropo_corr(filter(1, :)), x)');
iono_corr_gim = cat(3, iono_corr_gim, interp1(tmpIono_corr_gim(filter(1, :)), x)');
sea_state_bias = cat(3, sea_state_bias, interp1(tmpSea_state_bias(filter(1, :)), x)');

range = cat(3, range, interp1(tmpRange(filter), x)');

mss = cat(3, mss, interp1(tmpMss(filter(1,:)), x)');
ssha = cat(3, ssha, interp1(tmpSsha(filter(1,:)), x)');

solidEarthTideHeight = cat(3, solidEarthTideHeight, interp1(tmpSolidEarthTideHeight(filter(1,:)), x)');
oceanTide = cat(3, oceanTide, interp1(tmpOceanTide(filter(1,:)), x)');
poleTide = cat(3, poleTide, interp1(tmpPoleTide(filter(1,:)), x)');
invBarCorr = cat(3, invBarCorr, interp1(tmpInvBarCorr(filter(1,:)), x)');
HF = cat(3, HF, interp1(tmpHF(filter(1,:)), x)');

N = cat(1, N, tmpN);

n = n + 1;

% Remove tmp variables from workspace
clear -regexp ^tmp

%% Calculate height
% Hardware variables
C_ntp = 51;
B_spc = 0.31;

% Init calculated variables
C_rtrk_ocog = zeros(N, 1, n);
C_rtrk_pp_cog = zeros(N, 1, n);
M_ocog = zeros(N, 1, n);
W_ocog = zeros(N, 1, n);
M_pp_cog = zeros(N, 1, n);
W_pp_cog = zeros(N, 1, n);
mp = zeros(N, 1, n);
pP = zeros(N, 1, n);
wPower = zeros(128, N);

pStart = zeros(N, 1, n);
pStop = zeros(N, 1, n);
for i = 1:n
    for j = 1:N
        [C_rtrk_ocog(j,1,i), ~,~, M_ocog(j,1,i), W_ocog(j,1,i)] = waveformAnalysis(wave(:,j,i), 'OCOG');
        [C_rtrk_pp_cog(j,1,i), pStart(j,1,i), pStop(j,1,i), M_pp_cog(j,1,i), W_pp_cog(j,1,i)] = waveformAnalysis(wave(:,j,i), 'PP_COG');
        mp(j,1,i) = maxPower(wave(:,j,i), agc(j,1,i));
        pP(j,1,i) = pulsePeakiness(wave(:,j,i), 128);
        wPower(:,j) = wavePower(wave(:,j,i), agc(j,1,i));
    end
end

%% Correction application
% Retracked height calculated, both OCOG and PPCOG
epoch_ocog = (C_ntp - C_rtrk_ocog) * B_spc;
epoch_pp_cog = (C_ntp - C_rtrk_pp_cog) * B_spc;

% Tracker corrections
% Missing system bias
tracker_corr = modeled_instr_corr + doppler_corr;

% Corrected Altimerter Range
sea_state_bias(isnan(sea_state_bias)) = 0; % Sets NAN values to 0
alt_corr = model_dry_tropo_corr + rad_wet_tropo_corr + iono_corr_gim + ...
    sea_state_bias + tracker_corr;
correctedRange = tracker + alt_corr;

% Sea Surface Hegiht
ssh_ocog = alt - correctedRange + epoch_ocog;
ssh_pp_cog = alt - correctedRange + epoch_pp_cog;

% Sea height anomaly
sla_corr = mss + solidEarthTideHeight + oceanTide + poleTide + invBarCorr + HF;
sla_ocog = ssh_ocog - sla_corr;
sla_pp_cog = ssh_pp_cog - sla_corr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Our SLA OCOG vs product SLA
figure;
subplot(3, 1, 1);
plot(sla_ocog);
title('OCOG');

subplot(3, 1, 2);
plot(ssha);
title('Product ssha');

subplot(3, 1, 3);
plot((sla_ocog - ssha));
title('Difference');


%% Our SLA PPCOG vs product SLA
figure;
subplot(3, 1, 1);
plot(sla_pp_cog);
title('PPCOG');

subplot(3, 1, 2);
plot(ssha);
title('Product ssha');

subplot(3, 1, 3);
plot((sla_pp_cog - ssha));
title('Difference');

%% Waveforms OCOG
i = 1110;
figure;
subplot(2,1,1);
hold on;
plot(wave(:,i));
line([C_rtrk_ocog(i), C_rtrk_ocog(i)], get(gca,'ylim'), 'color', 'r', 'linestyle', '--');
line([C_ntp, C_ntp], get(gca,'ylim'), 'color', 'g', 'linestyle', '--');
rectangle('position', [C_rtrk_ocog(i) 0 W_ocog(i) M_ocog(i)]);
title('OCOG', 'fontSize', 18);
legend('wave', 'retracked point', 'reference bin')

% Waveforms PPCOG
subplot(2,1,2);
hold on;
plot(wave(:,i));
plot(pStart(i):pStop(i), wave(pStart(i):pStop(i),i), 'color', 'r');
line([C_rtrk_pp_cog(i), C_rtrk_pp_cog(i)], get(gca,'ylim'), 'color', 'r', 'linestyle', '--');
line([C_ntp, C_ntp], get(gca,'ylim'), 'color', 'g', 'linestyle', '--');
rectangle('position', [C_rtrk_pp_cog(i) 0 W_pp_cog(i) M_pp_cog(i)]);
title('PP COG (Threshold retracker)', 'fontSize', 18);
legend('wave', 'primary peak', 'retracked point', 'reference bin')

%%
% Waveforms PPCOG
figure;
hold on;
plot(wave(:,i));
plot(pStart(i):pStop(i), wave(pStart(i):pStop(i),i), 'color', 'r');
line([C_rtrk_pp_cog(i), C_rtrk_pp_cog(i)], get(gca,'ylim'), 'color', 'r', 'linestyle', '--');
line([C_ntp, C_ntp], get(gca,'ylim'), 'color', 'g', 'linestyle', '--');
rectangle('position', [C_rtrk_pp_cog(i) 0 W_pp_cog(i) M_pp_cog(i)], 'EdgeColor', [.25, .25, .25]);
plot(nan, nan, 'color', [.25, .25, .25])
title('PP COG (Threshold retracker)', 'fontSize', 18);
xlabel('Bin #');
ylabel('Count');
legend('wave', 'primary peak', 'retracked point', 'reference bin', 'OCOG rectangle');
fnam = sprintf('figures/ppcog_%s', cycleName);
hgexport(gcf, fnam, style);

%%
figure;
colormap(jet)
imagesc(wave())
set(gca, 'yDir', 'normal');
colorbar;
xlabel('Waveform #');
ylabel('Gate bin');
title('Image of track waveforms', 'fontSize', 18);
hgexport(gcf, strcat('figures\full_spectrum_' ,cycleName), style);

%%
figure;
xq = 1:N;
hold on;
plot(sla_pp_cog);
plot(xq(pP >= 30), sla_pp_cog(pP >= 30), 'r.');
xlim([0 N]);

%%
figure;
title('Pulse Peaknis, scatter');
m_scatter(lon, lat, 10, mp, 'filled');
m_gshhs('lc', 'color', 'k');
m_grid;

%%
i=1;
figure;
hold on;
plot(wave(:,i));
line([C_ntp, C_ntp], get(gca, 'ylim'), 'color', 'g', 'linestyle', '--');
title('Wave', 'fontSize', 18);
xlabel('Bin #');
ylabel('Counts');
set(gcf,'units','points','position',[50,250,595,420]);
% fnam = sprintf('figures\\single_wave_cycle_%s_N%f_E%f', cycleName, lat(i), lon(i));
% fnam = strrep(fnam, '.', '-');
% hgexport(gcf, fnam, style);

%%
if strcmp(input('Make movie (y/n)\n', 's'), 'y')
fig = figure('units', 'pixels', 'position', [400 100 1024 768]);
f = struct('cdata', cell(1, 100), 'colormap', cell(1, 100));
for i = 1:length(wave(1,:))
    subplot(2,1,1);
    plot(wave(:,i));
    line([C_rtrk_ocog(i), C_rtrk_ocog(i)], get(gca,'ylim'), 'color', 'r', 'linestyle', '--');
    line([C_ntp, C_ntp], get(gca,'ylim'), 'color', 'g', 'linestyle', '--');
    
    % Waveforms PPCOG
    subplot(2,1,2);
    plot(wave(:,i));
    line([C_rtrk_pp_cog(i), C_rtrk_pp_cog(i)], get(gca,'ylim'), 'color', 'r', 'linestyle', '--');
    line([C_ntp, C_ntp], get(gca,'ylim'), 'color', 'g', 'linestyle', '--');
    title(num2str(i));
    f(i) = getframe(fig);
end
end
