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
track = [[], []];
alt = [[], []];

modeled_instr_corr = [[], []];
doppler_corr = [[], []];

model_dry_tropo_corr = [[]  , []];
model_wet_tropo_corr = [[], []];
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
for cycle = 32
    cycleName = sprintf('cycle_%03d', cycle);
    cycleFile = strcat(cycleName, '.mat');
    
    if exist(cycleFile, 'file') == 0;
        % All data files
        cycleFilePath = fullfile(altikaFiles, cycleName);
        lis = dir(cycleFilePath);
        lis(1:2) = [];
        for j = 1:1002
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
                fprintf('%d skipped\n', j);
                continue;
            end
            
            % Tmp import of data variables
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
            track = [track; tmpTrack(filter)];
            alt = [alt; tmpAlt(filter)];
            
            modeled_instr_corr = [modeled_instr_corr; tmpModeled_intr_corr(filter(1, :))];
            doppler_corr = [doppler_corr; tmpDoppler_corr(filter(1,:))];
            
            model_dry_tropo_corr = [model_dry_tropo_corr; tmpModel_dry_tropo_corr(filter(1, :))];
            model_wet_tropo_corr = [model_wet_tropo_corr; tmpModel_wet_tropo_corr(filter(1, :))];
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
            
            disp(j)
        end
        % Save variables
        save(cycleFile, 'lon', 'lat', 'wave', 'agc', ...
            'track', 'alt', 'modeled_instr_corr', 'doppler_corr', ...
            'model_dry_tropo_corr', 'model_wet_tropo_corr', ...
            'iono_corr_gim', 'sea_state_bias', 'range', 'mss', 'ssha',...
            'solidEarthTideHeight', 'oceanTide', 'poleTide', ...
            'invBarCorr', 'HF');
    else
        % Load variables
        fprintf('Existing file %s has been loaded', cycleFile);
        load(cycleFile);
    end
end


%% Calculate height
% Hardware variables
C_ntp = 64;
B_spc = 0.31;

% Number of data points
N = length(wave(1,:));

% Init calculated variables
C_rtrk_ocog = zeros(N, 1);
C_rtrk_pp_cog = zeros(N, 1);
mp = zeros(N,1);
pP = zeros(N,1);

for i = 1:N
    C_rtrk_ocog(i) = waveformAnalysis(wave(:,i), 'OCOG');
    C_rtrk_pp_cog(i) = waveformAnalysis(wave(:,i), 'PP_COG');
    mp(i) = maxPower(wave(:,i), agc(i));
    pP(i) = pulsePeakness(wave(:,i), 128);
end

%% Classification
type = zeros(N, 1);
type(pP > 30) = 4;

type2 = zeros(N, 1);
% type2(40 <= pP < 50) = 2;
% type2(50 <= pP < 70) = 3;
type2(70 <= pP) = 4;

%% Correction application
x = linspace(1,N, length(modeled_instr_corr));

% Range corrections
track_corr = modeled_instr_corr + doppler_corr;
track_corr_40hz = interp1(x, track_corr, 1:N)';

% Altitude correction
sea_state_bias(isnan(sea_state_bias)) = 0; % Sets NAN values to 0
alt_corr = model_dry_tropo_corr + model_wet_tropo_corr + iono_corr_gim + sea_state_bias;
alt_corr_40hz = interp1(x, alt_corr, 1:N)';

% Sea height anomaly corrections
sla_corr = mss - solidEarthTideHeight - oceanTide - poleTide - invBarCorr - HF;
sla_corr_40hz = interp1(x, sla_corr, 1:N)';

% Retracked height calculated, both OCOG and PPCOG
epoch_ocog = (C_ntp - C_rtrk_ocog) * B_spc;
epoch_pp_cog = (C_ntp - C_rtrk_pp_cog) * B_spc;

% Sea surface height
ssh_ocog = alt + alt_corr_40hz - (track + track_corr_40hz) + epoch_ocog;
ssh_pp_cog = alt + alt_corr_40hz - (track + track_corr_40hz) + epoch_pp_cog;

% Sea surface height anamoly
sla = ssh_ocog - sla_corr_40hz;
sla_pp = ssh_pp_cog - sla_corr_40hz;

% Product sea surface height anamoly
ssha_40hz = interp1(x, ssha, 1:N)';

%% Plot
[Xq, Yq] = meshgrid(LON(1):0.1:LON(2), LAT(1):0.01:LAT(2));
Vq = griddata(lon, lat, sla_pp, Xq, Yq);

figure;
hold on
m_pcolor(Xq, Yq, Vq);
shading flat;
% m_scatter(lon, lat, 5, ssha);
m_gshhs('lc', 'color', 'k');
colorbar;
m_grid;

figure;
m_scatter(lon, lat, 10, sla_pp, 'filled');
m_gshhs('lc', 'color', 'k');
colorbar;
m_grid;

%%
Vq = griddata(lon, lat, type, Xq, Yq);

figure;
hold on
m_pcolor(Xq, Yq, Vq);
shading flat;
% m_scatter(lon, lat, 5, ssha);
m_gshhs('lc', 'color', 'k');
m_grid;

figure;
m_scatter(lon, lat, 10, type, 'filled');
m_gshhs('lc', 'color', 'k');
m_grid;