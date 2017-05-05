function result = ncImport(ncid, var)
varid = netcdf.inqVarID(ncid, var);
data = netcdf.getVar(ncid, varid, 'double');
fillValue = netcdf.getAtt(ncid, varid, '_FillValue');
data(data == fillValue) = NaN;
scale = 1;
offset = 0;
try
    scale = netcdf.getAtt(ncid, varid, 'scale_factor');
catch
end
try
    offset = netcdf.getAtt(ncid, varid, 'add_offset');
catch
end
result = data * scale + offset;
end