function [ velocities ] = velocity( fileName, LON, LAT )
%VELOCITY Extracts velocities from desired file
%   velocities = velocity(FILENAME, LON, LAT) extracts and parses
%   velocities from FILENAME and filters data within LON and LAT. LON and
%   LAT are arrays with start and stop value for area of interest.
% 
%   Example:
%       velocities = velocity('20160303.n.S1Adrift.vector', [-10, 10], [79,
%       82]);
%
%   See also GRIDVELOCITIES

%% Init. velocity struct
velocities = struct('lon', [], 'lat', [], 'x', [], 'y', [], 'magn', [], 'brng' ,[]);

%% Open file
fileID = fopen(fileName);
inputCell = textscan(fileID, '%f %f %f %f %s', 'Headerlines', 1);
fclose(fileID);

%% Make data useable
input = cell2mat(inputCell(1:4));

%% Extract data
skip = isnan(input(:,4));           % Skip human-readable line
start = input(skip,3) == 1;         % Start positions
stop = input(skip,3) == 3;          % Stop position

lon = input(skip,2);
lat = input(skip,1);

%% Filter for earth
filLon = (LON(1) - 1) < lon & lon < (LON(2) + 1) & stop;
filLat = (LAT(1) - 1) < lat & lat < (LAT(2) + 1) & stop;

filter = filLat & filLon;
filter = filter | circshift(filter, -1);

% Remove irrelevant data
lon(~filter) = [];
lat(~filter) = [];
start(~filter) = [];
stop(~filter) = [];

%% Populate velocities
velocities.lon = lon(stop);
velocities.lat = lat(stop);

velocities.x = coor2dist(lon(start), lat(start), lon(stop), lat(start));
velocities.y = coor2dist(lon(stop), lat(start), lon(stop), lat(stop));

% Assign direction 
velocities.x(lon(start) > lon(stop)) = velocities.x(lon(start) > lon(stop)) * -1;
velocities.y(lat(start) > lat(stop)) = velocities.y(lat(start) > lat(stop)) * -1;

velocities.brng = coor2brng(lon(start), lat(start), lon(stop), lat(stop));
velocities.magn = coor2dist(lon(start), lat(start), lon(stop), lat(stop));
end