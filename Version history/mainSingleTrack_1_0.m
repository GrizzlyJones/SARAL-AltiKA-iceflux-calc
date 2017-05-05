%% Init
clc; close all; clear;

%% File mangagement
addpath(fullfile(matlabroot, 'toolbox', 'matlab', 'm_map')); % m_maps
addpath(fullfile(pwd,'scripts'));                            % used scripts
altikaFiles = fullfile(pwd,'ALTIKA');                        % data

%% Fram Strait
LON = [-10, 10];
LAT = [76, 82];

% Settings for map projection
m_proj('albers equal-area', 'long', LON, 'lat', LAT, 'rectbox', 'off');

% Init. data arrays
lon = [[], []];
lat = [[], []];
wave = [[], []];
agc = [[], []];
tracker = [[], []];
alt = [[], []];

modeled_instr_corr = [[], []];
doppler_corr = [[], []];

model_dry_tropo_corr = [[]  , []];
rad_wet_tropo_corr = [[], []];
iono_corr_gim = [[], []];
sea_state_bias = [[], []];

range = [[], []];

mss = [[], []];
ssha = [[], []];

solidEarthTideHeight = [[], []];
oceanTide = [[], []];
poleTide = [[], []];
invBarCorr = [[], []];
HF = [[], []];

%% Load data
cycle = 32;
j = 18;

cycleName = sprintf('cycle_%03d', cycle);
cycleFile = strcat(cycleName, '.mat');
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
if ~(any(any(filter)))
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

% Filtration of found data, saving for later use
lon = [lon; tmpLon(filter)];
lat = [lat; tmpLat(filter)];
wave = [wave, tmpWave(:,filter)];
agc = [agc; tmpAGC(filter)];
tracker = [tracker; tmpTracker(filter)];
alt = [alt; tmpAlt(filter)];

modeled_instr_corr = [modeled_instr_corr; tmpModeled_intr_corr(filter(1, :))];
doppler_corr = [doppler_corr; tmpDoppler_corr(filter(1,:))];

model_dry_tropo_corr = [model_dry_tropo_corr; tmpModel_dry_tropo_corr(filter(1, :))];
rad_wet_tropo_corr = [rad_wet_tropo_corr; tmpRad_wet_tropo_corr(filter(1, :))];
iono_corr_gim = [iono_corr_gim; tmpIono_corr_gim(filter(1, :))];
sea_state_bias = [sea_state_bias; tmpSea_state_bias(filter(1, :))];

range = [range; tmpRange(filter)];

mss = [mss; tmpMss(filter(1,:))];
ssha = [ssha; tmpSsha(filter(1,:))];

solidEarthTideHeight = [solidEarthTideHeight; tmpSolidEarthTideHeight(filter(1,:))];
oceanTide = [oceanTide; tmpOceanTide(filter(1,:))];
poleTide = [poleTide; tmpPoleTide(filter(1,:))];
invBarCorr = [invBarCorr; tmpInvBarCorr(filter(1,:))];
HF = [HF; tmpHF(filter(1,:))];

%% Calculate height
% Hardware variables
C_ntp = 51;
B_spc = 0.31;

% Number of data points
N = length(wave(1,:));

% Init calculated variables
C_rtrk_ocog = zeros(N, 1);
C_rtrk_pp_cog = zeros(N, 1);
mp = zeros(N,1);
pP = zeros(N,1);

pStart = zeros(N,1);
pStop = zeros(N,1);

for i = 1:N
    C_rtrk_ocog(i) = waveformAnalysis(wave(:,i), 'OCOG');
    [C_rtrk_pp_cog(i), pStart(i), pStop(i)] = waveformAnalysis(wave(:,i), 'PP_COG');
    mp(i) = maxPower(wave(:,i), agc(i));
    pP(i) = pulsePeakness(wave(:,i), 128);
end

%% Correction application
x = linspace(1,N, length(modeled_instr_corr));

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
alt_corr_40hz = interp1(x, alt_corr, 1:N)';
correctedRange = tracker + alt_corr_40hz;

% Sea Surface Hegiht
ssh_ocog = alt - correctedRange + epoch_ocog;
ssh_pp_cog = alt - correctedRange + epoch_pp_cog;

% Sea height anomaly
sla_corr = mss + solidEarthTideHeight + oceanTide + poleTide + invBarCorr + HF;
sla_corr_40hz = interp1(x, sla_corr, 1:N)';
sla_ocog = ssh_ocog - sla_corr_40hz;
sla_pp_cog = ssh_pp_cog - sla_corr_40hz;

% Product sea surface height anamoly
ssha_40hz = interp1(x, ssha, 1:N)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Our SLA OCOG vs product SLA
figure;
subplot(3, 1, 1);
plot(sla_ocog);
title('OCOG');

subplot(3, 1, 2);
plot(ssha_40hz);

subplot(3, 1, 3);
plot(sla_ocog - ssha_40hz);


%% Our SLA PPCOG vs product SLA
figure;
subplot(3, 1, 1);
plot(sla_pp_cog);
title('PPCOG');

subplot(3, 1, 2);
plot(ssha_40hz);

subplot(3, 1, 3);
plot(sla_pp_cog - ssha_40hz);

%% Waveforms OCOG
figure;
subplot(2,1,1);
hold on;
plot(wave(:,1));
line([C_rtrk_ocog(1), C_rtrk_ocog(1)], get(gca,'ylim'), 'color', 'r');
line([C_ntp, C_ntp], get(gca,'ylim'), 'color', 'g');

% Waveforms PPCOG
subplot(2,1,2);
hold on;
plot(wave(:,1));
plot(pStart:pStop, wave(pStart:pStop,1), 'color', 'r');
line([C_rtrk_pp_cog(1), C_rtrk_pp_cog(1)], get(gca,'ylim'), 'color', 'r');
line([C_ntp, C_ntp], get(gca,'ylim'), 'color', 'g');

%%
fig = figure('units', 'pixels', 'position', [400 100 1024 768]);
for i = 1:N
    subplot(2,1,1);
    plot(wave(:,i));
    line([C_rtrk_ocog(i), C_rtrk_ocog(i)], get(gca,'ylim'), 'color', 'r');
    line([C_ntp, C_ntp], get(gca,'ylim'), 'color', 'g');
    
    % Waveforms PPCOG
    subplot(2,1,2);
    plot(wave(:,i));
    line([C_rtrk_pp_cog(i), C_rtrk_pp_cog(i)], get(gca,'ylim'), 'color', 'r');
    line([C_ntp, C_ntp], get(gca,'ylim'), 'color', 'g');
    title(num2str(i));
    f(i) = getframe(fig);
end
