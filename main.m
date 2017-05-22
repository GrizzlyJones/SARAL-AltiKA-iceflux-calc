%% Init
clc; close all; clear;

%% File mangagement
addpath(fullfile(matlabroot, 'toolbox', 'matlab', 'm_map')); % m_maps
addpath(genpath(fullfile(pwd,'scripts')));                   % used scripts
altikaFiles = 'D:\Altika';                                   % data

%% Fram Strait
LON = [-10, 10];
LAT = [79, 82];

% Settings for map projection
m_proj('albers equal-area', 'long', LON, 'lat', LAT, 'rectbox', 'off');

% Define fluxgate
fluxgate = initFluxgate([-8.2, 8.9], [81.4, 80], 1000);

% Init. data arrays
data = struct('lon', [], 'lat', [], ...
    'wave', [], 'agc', [], 'tracker', [], 'alt', [], ...
    'modeled_instr_corr', [], 'doppler_corr', [], ...
    'model_dry_tropo_corr', [], 'rad_wet_tropo_corr', [], 'iono_corr_gim', [], 'sea_state_bias', [], ...
    'range', [], 'mss', [], 'ssha', [], ...
    'solidEarthTideHeight', [], 'oceanTide', [], 'poleTide', [], 'invBarCorr', [], 'HF', []);
dataNames = fieldnames(data);
for name = dataNames
    data.(name{1}) = [];
end
dataNames(1:2)=[];
names40hz = {'waveforms_40hz', 'agc_40hz', 'tracker_40hz', 'alt_40hz'};
corrNames = {'modeled_instr_corr_range', 'doppler_corr', ...
    'model_dry_tropo_corr', 'rad_wet_tropo_corr', 'iono_corr_gim', ...
    'sea_state_bias', 'range_40hz', 'mean_sea_surface', 'ssha', ...
    'solid_earth_tide', 'ocean_tide_sol2', 'pole_tide', ...
    'inv_bar_corr', 'hf_fluctuations_corr'};

%% Load data
for cycle = 32
    cycleName = sprintf('cycle_%03d', cycle);
    cycleFile = fullfile(pwd,'data', strcat(cycleName, '.mat'));
    
    if exist(cycleFile, 'file') == 0
        disp('No file found, creating new.');
        % All data files
        cycleFilePath = fullfile(altikaFiles, cycleName);
        lis = dir(cycleFilePath);
        lis(1:2) = [];
        %         for j = 1:length(lis)
        for j = 1:length(lis)
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
                % fprintf('%d skipped\n', j);
                continue;
            end
            
            % Filtration of found data, saving for later use
            data.lon = vertcat(data.lon, tmpLon(filter));
            data.lat = vertcat(data.lat, tmpLat(filter));
            
            for k = 1:length(names40hz)
                disp(dataNames{k})
               data.(dataNames{k}) = load40hz(data.(dataNames{k}), filePath, names40hz{k}, filter);
            end
            
            for k = 1:length(corrNames)
                disp(dataNames{k+4})
               data.(dataNames{k+4}) = loadCorr(data.(dataNames{k+4}), filePath, corrNames{k}, filter); 
            end
            
            % Remove tmp variables from workspace
            clear -regexp ^tmp
                        
            % disp(j)
        end
        
        % Save variables
        save(cycleFile, 'data');
    else
        % Load variables
        fprintf('Existing file %s has been loaded\n', cycleFile);
        load(cycleFile);
    end
    
    % Load velocities
    velocities = velocity('velocity\20160303.n.S1Adrift.vector', LON, LAT);
end

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
W = zeros(N,1);

for i = 1:N
    C_rtrk_ocog(i) = waveformAnalysis(wave(:,i), 'OCOG');
    [C_rtrk_pp_cog(i), ~, ~, ~, W(i)] = waveformAnalysis(wave(:,i), 'PP_COG');
    mp(i) = maxPower(wave(:,i), agc(i));
    pP(i) = pulsePeakness(wave(:,i), 128);
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

%% Grid iterpolation
gridData = struct('Xq', [], 'Yq', [], 'sla', [], 'ssha', [], ...
                  'pP', [], 'mP', [], 'W', []);
gridData.sla = struct('ocog', [], 'pp_cog', []);

[gridData.Xq, gridData.Yq] = meshgrid(LON(1):0.01:LON(2), LAT(1):0.001:LAT(2));
gridData.sla.pp_cog = griddata(lon, lat, sla_pp_cog, Xq, Yq);
gridData.ssha = griddata(lon, lat, ssha, Xq, Yq);
gridData.pP = griddata(lon, lat, pP, Xq, Yq);
gridData.mP = griddata(lon, lat, mp, Xq, Yq);
gridData.W = griddata(lon, lat, W, Xq, Yq);

% Classification
pP_class = zeros(size(Xq));
pP_class(gridData.pP >= 30 & gridData.W < 2) = 4;

mP_class = zeros(size(Xq));
mP_class(gridData.mP >= 70) = 4;

%% Mask
mask = ~isnan(gridData.ssha);

%% Plot
% Primary Peak COG Map
figure;
hold on
title('Retracked SLA');
m_pcolor(Xq, Yq, gridData.sla.pp_cog);
shading flat;
colorbar;
m_grid;

% Primary Peak COG Map, scatter
figure;
title('Retracked, scatter');
m_scatter(lon, lat, 10, sla_pp_cog, 'filled');
colorbar;
m_grid;

%% Pulse Peakniss and Max Power
figure;
subplot(2,1,1);
title('Pulse Peakniss');
hold on
m_pcolor(Xq, Yq, pP_class);
shading flat;
m_grid;

subplot(2,1,2);
title('Max Power');
hold on
m_pcolor(Xq, Yq, mP_class);
shading flat;
m_grid;

%% Product given SLA
figure;
title('Product SLA');
hold on
m_pcolor(Xq, Yq, gridData.ssha);
shading flat;
colorbar;
m_grid;

%% Track grid vs interp
steps = 1000;
lonx = linspace(-8.2, 8.9, steps);
latx = linspace(81.4, 80, steps);
fluxgate_sla_pp_cog = interp2(Xq, Yq, gridData.sla.pp_cog, lonx, latx);
fluxgate_pP = interp2(Xq, Yq, gridData.pP, lonx, latx);
fluxgate_W = interp2(Xq, Yq, gridData.W, lonx, latx);

figure
hold on
m_pcolor(Xq, Yq, gridData.ssha);
m_quiver(velocities.lon, velocities.lat, velocities.x, velocities.y);
shading flat;
colorbar;
m_grid;
title('Ice drift');

%% Track grid vs interp
fluxgate = interpProfile(fluxgate, gridData, gridVelocity);

freeboard = freeboardAnalysis(fluxgate);
freeboard = thickness(freeboard, 'radar');
normVelocities = projVelocity(fluxgate);

[volFlow, flow] = calcVolFlow(fluxgate, freeboard, normVelocities);
fprintf('The volumetric flow is %.2f cubic kilometers per day (positive is northen flow)\n', volFlow*1e-9);

figure;
plot(fluxgate.profile.cumStep, flow);

plotFluxgate(gridData, fluxgate, freeboard);