function [ data ] = loadData( array, filePath, variable, filter )
%LOADDATA Summary of this function goes here
%   Detailed explanation goes here


names40hz = {'waveforms_40hz', 'agc_40hz', 'tracker_40hz', 'alt_40hz'};
corrNames = {'modeled_instr_corr_range', 'doppler_corr', ...
    'model_dry_tropo_corr', 'rad_wet_tropo_corr', 'iono_corr_gim', ...
    'sea_state_bias', 'range_40hz', 'mean_sea_surface', 'ssha', ...
    'solid_earth_tide', 'ocean_tide_sol2', 'pole_tide', ...
    'inv_bar_corr', 'hf_fluctuations_corr'};

tmpData = ncread(filePath, variable);

if ~isempty(strfind(variable, 'waveforms'))
    data = vertcat(array, tmpData(:, filter));
elseif ~isempty(strfin(variable, '_40hz'))
    data = vertcat(array, tmpData(filter));
else
    N = sum(filter(:));
    n = sum(filter(1,:));

    x = linspace(1, n, N);

    data = vertcat(array, interp1(tmpData(filter(1, :)), x)');
end

end

