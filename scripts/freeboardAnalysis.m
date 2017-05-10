function iceSheets = freeboardAnalysis(fluxgate)
iceSheets = struct('lon', [], 'lat', [], 'freeboard', [], 'index', []);
iceSheets.index = struct('start', [], 'stop', []);

tmp = struct('heights', [], 'leads', []);

lon = fluxgate.lon;
lat = fluxgate.lat;
heights = fluxgate.sla;
leads = fluxgate.class;

if length(heights) ~= length(leads)
    error('Error:\nDimensions don\'' match');
end

N = length(heights);

initial = true;
newFreeboard = false;
index = 1;
indexStop = nan;

for i = 1:N
    if leads(i) == 4 && ~newFreeboard
        tmpStop = indexStop;
        indexStop = i - 1;
        newFreeboard = true;
        
        if ~initial
            if any(isnan(heights(indexStart:indexStop))) || ((indexStop - indexStart) < 3)
                continue;
            end
            iceSheets(index).lon = lon(indexStart:indexStop);
            iceSheets(index).lat = lat(indexStart:indexStop);
            iceSheets(index).index.start = indexStart;
            iceSheets(index).index.stop = indexStop;
            
            tmp(index).heights = heights(indexStart:indexStop);
            tmp(index).leads = heights(tmpStop:indexStart);
            
            index = index + 1;
        end
    elseif newFreeboard
        indexStart = i;
        newFreeboard = false;
        initial = false;
    end
end
tmp(index).leads = heights(tmpStop:indexStart);

for i = 1:length(iceSheets)
    avgLead = (mean(tmp(i).leads) + mean(tmp(i+1).leads)) / 2;
    iceSheets(i).freeboard = tmp(i).heights - avgLead;
end
end