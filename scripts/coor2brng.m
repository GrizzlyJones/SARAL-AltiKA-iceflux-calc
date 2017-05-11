function [ brng ] = coor2brng( lon1, lat1, lon2, lat2 )
%COOR2BRNG The bearing from coordinat 1 to coordinat 2
%   brng = coor2brng(LON1, LAT1, LON2, LAT2) calculates the bearing,
%   heading from LON1, LAT1, to LON2, LAT2.
%
%   Source: www.movable-type.co.uk/acripta/latlong.html
%           by Chris Veness

phi1 = deg2rad(lat1);
phi2 = deg2rad(lat2);

DeltaLambda = deg2rad(lon2 - lon1);

y = sin(DeltaLambda) .* cos(phi2);
x = cos(phi1) .* sin(phi2) - sin(phi1) .* cos(phi2) .* cos(DeltaLambda);

brng = rad2deg(atan2(y, x));

end

