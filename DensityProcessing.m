% Author: Riley Culberg
% Date: 12/7/2020
%
% This scripts take inverted radar reflectivities, extracts the 2012 melt
% layer, and averages over 1km bins to create the final density dataset
% shown in Figure 2.
%
% Data dependencies:
% Inverted radar reflectivities - /DerivedData/Demo_Density.mat
%   -> This input data file can also be created by running the top-level
%   script called DensityInversion.m 
% Mean climate data - /RefData/Greenland_5km_v1.1.nc
% ------------------------------------------------------------------------

clear;

addpath(genpath('./ReferenceFunctions/'));

% Calculate the approximate expected depth of the 2012 summer surface using
% the steady state Herron & Langway model - to be used later for rough
% localization of the melt layer in the radar data

rho_ice = 0.917;          % density of solid ice in g/cm^3
R = 8.314;                % gas constant
rho1 = 0.33:0.001:0.55;   % density axis

% Load mean annual accumulation rates and surface temperatures
accum = ncread('./RefData/Greenland_5km_v1.1.nc', 'presprcp');
temp = ncread('./RefData/Greenland_5km_v1.1.nc', 'surftemp');

% Convert accumulation from mWe to mIe
accum = accum*(997/917);

depth = zeros(size(accum));
for m = 1:size(accum,1)
    for n = 1:size(accum,2)
        
        A = accum(m,n);
        T = 273.15 + temp(m,n);
        
        % Herron & Langway steady state solution for depth and age over the
        % firn in the 1st compaction regime
        k0 = 11*exp(-10160/(R*T));
        h1 = (1./(rho_ice.*k0)).*(log(rho1./(rho_ice-rho1)) - log(rho1(1)./(rho_ice - rho1(1))));
        t1 = (1/(k0*A)).*log((rho_ice - rho1(1))./(rho_ice - rho1));
        
        % Calculate the depth of the 2012 horizon in 2017 (~4.66 years
        % after initial deposition)
        [~, ind] = min(abs(t1 - 4.66));
        depth(m,n) = h1(ind);
    end
end

% Load the inverted radar reflectivities
load('./DemoData/Demo_Density.mat');

% Create a single combined data structure with results from every radargram
% in the flight transect
data_min = [];
data_max = [];
lat = [];
lon = [];
for k = 1:length(density_result)
    data_min = [data_min density_result{k}.inversion_min];
    data_max = [data_max density_result{k}.inversion_max];
    lat = [lat density_result{k}.lat]; 
    lon = [lon density_result{k}.lon]; 
end

% Calculate the approximate depth of the 2012 melt layer of this flight
% transect
lat_ref = ncread('./RefData/Greenland_5km_v1.1.nc', 'lat');
lon_ref = ncread('./RefData/Greenland_5km_v1.1.nc', 'lon');

line_depth = zeros(size(lat));
for m = 1:length(line_depth)
        opt = abs(lat_ref - lat(m)).^2 + abs(lon_ref - lon(m)).^2;
        [value, index] = min(opt(:));
        [row, col] = ind2sub(size(opt), index);
        line_depth(m) = depth(row,col);
end

% Find the minimum and maximum consistent density for each trace with +/- 8
% samples of the estimated layer depth
min_density = zeros(1, size(data_min, 2));
max_density = zeros(1, size(data_min, 2));
for k = 1:length(min_density)
    
    % Establish layer extraction vertical extent
    start = round((line_depth(k)/0.3) - 8);
    stop = round((line_depth(k)/0.3) + 8);
    if start < 1
        start = 1;
    end
    if stop > 50
        stop = 50;
    end
    
    % Pull min and max values, ignoring traces set to NaN where data was
    % discarded due to high roll or 1000 where no inversion was run due to
    % no layer detection
    seg = data_min(start:stop,k);
    seg2 = data_max(start:stop,k);
    x = find(~isnan(seg));
    if seg(1) ~= 1000 && ~isempty(seg(x))
        min_density(k) = max(seg(x));
        max_density(k) = max(seg2(x));
    elseif seg(1) == 1000
        min_density(k) = 1000;
        max_density(k) = 1000;
    else
        min_density(k) = NaN;
        max_density(k) = NaN;
    end
end

% Average the layer density over 1km bins
[at, ~, ~, ~] = geodetic_to_along_track(lat, lon, []);
x = 0:1000:at(end);

mean_lat = zeros(size(x));
mean_lon = zeros(size(x));
mean_min_density = NaN*ones(size(x));
mean_max_density = NaN*ones(size(x));
coverage = zeros(size(x));
for m = 1:length(x) - 1
    [~, start] = min(abs(at - x(m)));
    [~, stop] = min(abs(at - x(m + 1)));
    mean_lat(m) = mean(lat(start:stop));
    mean_lon(m) = mean(lon(start:stop));
    tmp1 = min_density(start:stop);
    tmp2 = max_density(start:stop);
    q = find(~isnan(tmp1));
    tmp1 = tmp1(q);
    tmp2 = tmp2(q);
    p = find(tmp1 == 1000);
    tmp1(p) = [];
    tmp2(p) = [];
    coverage(m) = (length(q) - length(p))/(stop - start + 1 - length(p));
    if ~isempty(tmp1) && stop > start + 1
        mean_min_density(m) = mean(tmp1);
        mean_max_density(m) = mean(tmp2);
    end
end

% Remove any traces where the density is NaN due to excluded data
ind = [];
for k = 1:length(mean_min_density)
    if isnan(mean_min_density(k))
        ind = [ind k];
    end
end

mean_min_density(ind) = [];
mean_max_density(ind) = [];
coverage(ind) = [];
mean_lat(ind) = [];
mean_lon(ind) = [];

% Write data to a text file
data = [mean_lat' mean_lon' mean_min_density' mean_max_density' coverage'];
writematrix(data, './DerivedData/Demo_Density_Coverage.txt');

