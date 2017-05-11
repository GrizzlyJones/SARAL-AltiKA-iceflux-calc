function [ lon2, lat2 ] = addDist2coor( lon1, lat1, brng, dist )
%ADDDIST2COOR Adds a distance in a certain bearing to a coordinat
%   [LON2, LAT2] = addDist2coord(LON1, LAT1, BRNG, DIST) adds distance,
%   DIST, to coordinate LON1, LAT1, in the heading of the bearing, BRNG.
%
%   Source: www.movable-type.co.uk/acripta/latlong.html
%           by Chris Veness

R = 6.371e6;

phi1 = deg2rad(lat1);
lambda1 = deg2rad(lon1);
theta = deg2rad(brng);
delta = dist / R;

phi2 = asin(sin(phi1) * cos(dela) + cos(phi1) * sin(delta) * cos(theta));
lambda2 = lambda1 + atan2(sin(theta) * sin(delta) * cos(phi1), cos(delta) - sin(phi1) * sin(phi2));

lon2 = rad2deg(lambda2);
lat2 = rad2deg(phi2);
end

