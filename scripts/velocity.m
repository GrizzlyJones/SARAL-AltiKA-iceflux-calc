function [ velocities ] = velocity( fileName, LON, LAT )
%VELOCITY extracts velocities from desired file
%   
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
filLon = (LON(1) - 1) < lon & lon < (LON(2) + 1) & stop;
lat = input(skip,1);
filLat = (LAT(1) - 1) < lat & lat < (LAT(2) + 1) & stop;

%% Filter for earth
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