function [ fluxgate ] = initFluxgate( lon, lat, steps )
%INITFLUXGATE Initiates a fluxgate with specified path and steps
%   fluxgate = initFluxgate(LON, LAT, STEPS, XQ, YQ, SLA_PPCOG_Q, SSHA_Q, PPQ, WQ, GRIDVELOCITY)
%   initiates a fluxgate from LON and LAT, with a detail of STEPS.
% 
%   Structure:
%   fluxgate
%       profile
%           lon         longitude
%           lat         latitude
%           brng        bearing
%           steps       steps between start and stop
%           stepSize    width of each step
%       data
%           sla         sea level anomaly
%           ssha        product sea leavel anomaly
%           magn        velocity magnitude
%           brng        velocity bearing
%           pP          pulse peakness
%           W           wave width
%           class       wave classification
%
%   See also FREEBOARDANALYSIS, INTERPPROFILE

%% Generate struct
fluxgate = struct('profile', [], 'data', []);
fluxgate.profile = struct('lon', [], 'lat', [], 'brng', [], ...
                          'steps', steps, 'stepSize', []);
fluxgate.data = struct('sla', [], 'ssha', [], 'magn', [], 'brng', [], ...
                       'pP', [], 'W', [], 'class', []);

%% Creating path for fluxgate
fluxgate.profile.lon = linspace(lon(1), lon(2), steps);
fluxgate.profile.lat = linspace(lat(1), lat(2), steps);
fluxgate.profile.brng = coor2brng(lon(1), lat(1), lon(2), lat(2));
fluxgate.profile.stepSize = coor2dist(lon(1), lat(1), lon(2), lat(2)) / steps;
end

