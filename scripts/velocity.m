function [t] = velocity(fileName)
% Indhentning af hastighederne fra www.seaice.dk/N/YYYY/MM/DD/

t = struct('lon', [], 'lat', [], 'x', [], 'y', [], 'magn', [], 'brng' ,[] , 'xy', [], 'XY', []);

fileID = fopen(fileName);
M = textscan(fileID, '%f %f %f %f %s', 'Headerlines', 1);
fclose(fileID);

A = cell2mat(M(1:4));
skip = isnan(A(:,4));
start = A(skip,3);

lon = A(skip,2);
lat = A(skip,1);

t.lon = lon(start == 3);
t.lat = lat(start == 3);

t.x = coor2dist(lon(start == 1), lat(start == 3), lon(start == 3), lat(start == 3));
t.y = coor2dist(lon(start == 1), lat(start == 1), lon(start == 1), lat(start == 3)); 

t.xy = sqrt(t.x.^2 + t.y.^2);
t.XY = coor2dist(lon(start == 1), lat(start == 1), lon(start == 3), lat(start == 3)); 


t.brng = coor2brng(lon(start == 1), lat(start == 1), lon(start == 3), lat(start == 3));

t.magn = coor2dist(lon(start == 1), lat(start == 1), lon(start == 3), lat(start == 3));

% B = cell2mat(M(4)); 
% C = celldisp(M(5));

end 

