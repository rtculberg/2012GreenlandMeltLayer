% Author: Riley Culberg
% Date: 12/7/2020
%
% This scripts calculates the absolute calibration coefficients for flight
% lines that intersect the North Greenland Traverse firn core sites (see
% Horhold, et al 2011 for discussion of these high resolution firn core
% measurements). As long as this flight line is also part of the leveling 
% network used in CrossSeasonLeveling.m and those offsets have been 
% normalized such that the offset at this tie transect is 0, then this 
% output value can be added to the coefficients produced by that script to 
% get the final calibration constants for all transects.
%
% Outputs:
%   S_final - absolute calibration coefficient in linear power units
%   (scalar double)
%
% Data dependencies:
%   High resolution firn core density profiles - ./RefData/density_268.mat 
%   include to demonstrate functionality
%   CReSIS L1B radargrams at the firn core sites
%       ./DemoData/Data_20120330_01_106.mat included to demonstrate the
%       functionality
%
% Function dependencies:
%   ./ReferenceFunctions/LayerSimulator_UWB2.m
% ------------------------------------------------------------------------
clear;

addpath(genpath('./ReferenceFunctions/'));

% Load the high resolution firn core measurements from B
load('./RefData/density_269.mat'); 

fc = 750e6;           % Radar center frequency (Hz)
BW = 300e6;           % Radar bandwidth (Hz)
tau = 2.048e-6;       % Radar pulse length (s)

% Convert density to permittivity with Kovacs relation
n = 1 + 0.845*density;
n = vertcat(1, n);
n_depth = midpoint_depth;
n_depth = vertcat(0, n_depth);
sampling = 0.3;       % Radar depth sampling
bin = 0;              % Set all to zero since not adding ice lenses
width = 0;
ice_density = 0;

% Simulate the expected firn reflection coefficients at this location
[layers, layer_depth] = LayerSimulator_UWB2(fc, BW, tau, n, n_depth, sampling, bin, width, ice_density);

% Load the radar data 
data1 = load('./DemoData/Data_20120330_01_106.mat');
data1.Data = 10*log10(data1.Data);

% Retrack the surface

c = 299792458;         % speed of light in a vacuum

start = 1528;           % Bound the indices where the surface is
stop = 1560;

% Find the max gradient within the surface bound
surf_ind = zeros(size(data1.Latitude));
for k = 1:length(surf_ind)
    max_del_y = 0;
    surface = 0;
    for m = start:1:stop
        del_y = data1.Data(m,k) - data1.Data(m-1,k);
        if del_y > max_del_y
            max_del_y = del_y;
            surface = m;
        end
    end
    surf_ind(k) = surface;
end

% Smooth surface along-track
for k = 2:length(surf_ind)
    if abs(surf_ind(k) - surf_ind(k-1)) > 3
        surf_ind(k) = surf_ind(k-1);
    end
end

% Find the first maximum after the smoothed max gradient and record the
% surface power and surface clearance
h = zeros(size(surf_ind));
for k = 1:length(surf_ind)
    surface = surf_ind(k);
    surf_power = data1.Data(surface,k);
    for m = surface+1:1:(surface + 4)
        if data1.Data(m,k) > surf_power
            surf_power = data1.Data(m,k);
            surface = m;
        end
        if data1.Data(m+1,k) < data1.Data(m,k)
            break;
        end
    end
    surf_ind(k) = surface;
    h(k) = 0.5*c*data1.Time(surface);
end

% Firn core location
lat = 77.2533;       
lon = -49.2167;

% Convert radar data to linear power units
data1.Data = 10.^(data1.Data/10);

% Find the trace nearest the firn core
[x1, y1] = ll2psn(lat, lon);
[x, y] = ll2psn(data1.Latitude, data1.Longitude);
[dist, ilat1] = min(abs(x - x1).^2 + abs(y- y1).^2);

%%%%%%%%%%%%%%%%%%% Calibrate the radar data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

num = 30;    % Number of surrounding traces to calibrate

surface = surf_ind(ilat1-num:ilat1+num);
h = h(ilat1-num:ilat1+num);

% True two-way travel time axis accounting for variable index of refraction
t = cumsum(2*(vertcat(0, diff(n_depth)).*n)/c);

S = zeros(size(surface));
for k = 1:length(surface)
    
    % Convert radar TWTT to depth
    t_true = data1.Time - data1.Time(surface(k));
    t_true = t_true(surface(k):end);
    
    depth1 = zeros(size(t_true));
    for m = 1:length(depth1)
        [~, ind] = min(abs(t - t_true(m)));
        depth1(m) = n_depth(ind);
    end
    
    % Set up geometric spreading loss correction
    depth1 = unique(depth1);
    geom = (2*(h(k) + depth1/1.6)).^2;
    
    % Geometrically correct radar received power and clip both radar and
    % simulation traces to depths between 15 and 80 meters to avoid any
    % influence from near-surface percolation
    x1 = 15:0.1:80;
    radar1 = interp1(depth1, data1.Data(surface(k):(surface(k)-1)+length(depth1),ilat1 - num - 1 + k).*geom, x1);
    sim1 = interp1(layer_depth, layers, x1);
    
    % Find the DC offset that minimizes the error between radar and
    % simulation in the mean square sense
    S(k) = radar1'\sim1';
end

S_final = mean(S);

save('./DerivedData/CalibrationConstant.mat', 'S_final');
