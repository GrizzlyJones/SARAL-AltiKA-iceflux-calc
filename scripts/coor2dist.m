function [ dist ] = coor2dist( lon1, lat1, lon2, lat2 )
%COOR2DIST Distance between two coordinates
%   dist = coor2dist(LON1, LAT1, LON2, LAT2) calculates the distance
%   between LON1, LAT1, and LON2, LAT2. The function uses the 'haversine'
%   formula to calculate the distances.
%   
%   Source: www.movable-type.co.uk/acripta/latlong.html
%           by Chris Veness

R = 6.371e6;

phi1 = deg2rad(lat1);
phi2 = deg2rad(lat2);

deltaPhi = deg2rad(lat2 - lat1);
deltaLambda = deg2rad(lon2 - lon1);

a = sin(deltaPhi/2) * sin(deltaPhi/2) + cos(phi1) * cos(phi2) * sin(deltaLambda/2) * sin(deltaPhi/2);
c = 2 * atan2(sqrt(a), sqrt(1-a));

dist = R * c;
end

