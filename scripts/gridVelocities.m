function [ gridVelocity ] = gridVelocities( Xq, Yq, velocities )
%GRIDVELOCITIES interpolates velocites to a grid
%   gridVelocity = gridVelocities(XQ, YQ, VELOCITIES) takes the bearing
%   (VELOCITIES.BRNG) and the magnintude (VELOCITIES.MANG) of the
%   velocities given in VELOCITIES and interpolates it to the given grid.
%
%   See also VELOCITY, INITFLUXGATE

warning('Implamentation of grid data will change in future versions');

gridVelocity = struct('brng', [], 'magn', []);

gridVelocity.brng = griddata(velocities.lon, velocities.lat, velocities.brng, Xq, Yq);
gridVelocity.magn = griddata(velocities.lon, velocities.lat, velocities.magn, Xq, Yq);
end

