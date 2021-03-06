% Author: Riley Culberg
% Date: 12/7/2020
%
% This scripts takes layer detections and inverts the received radar power 
% to estimate the minimum and maximum layer densities consistent with those
% observations. 
%
% Outputs:
% [Date]_Density.mat file for each flight. This is a 1xP cell array with one
% cell for each radargrams (segments) in the flight transect. Each cell
% contains the following data structures:
%    inversion_min - MxN matrix of doubles, where M is the fast time axis 
%    and N is the along-track axis. This contains the minimum density 
%    consistent with the radar observations. It takes a value of NaN if 
%    no layer was detected and 1000 if data was discarded due to excessive 
%    aircraft roll. 
%    inversion_max - MxN matrix of doubles, where M is the fast time axis 
%    and N is the along-track axis. This contains the maximum density 
%    consistent with the radar observations. It takes a value of NaN if 
%    no layer was detected and -1000 if data was discarded due to excessive 
%    aircraft roll. 
%    lat - N vector of latitude of each trace (WGS84)
%    lon - N vector of longitude for each trace (WGS84)
%    elev - N vector of surface elevation for each trace (WGS84)
%    base_density - MxN matrix of doubles, where M is the fast time axis 
%    and N is the along-track axis. This contains the background density 
%    assumed for each sample.
%
% Data dependencies:
% Background firn density - /RefData/BackgrounDensity.mat
%   - This reference data can also be generated with the
%      RunBackgroundDensity.m script
% Radar calibration constants - /DemoData/2017CalibrationConstants.mat
%   - These calibration constants were originally generated by running
%     Block 1 of the scripts. 
% Layer detections - /DemoData/Demo_Layer_Detections
% Precalculated radiometric references -
% /RefData/RadiometricReferenceCube.mat
%   - The radiometric reference cube was generated by running
%     RadiometricReference.mat
% 
% Function dependencies:
% ./ReferenceFunctions/C2xyz.m
% ------------------------------------------------------------------------

clear;

addpath(genpath('./ReferenceFunctions/'));

% Load Reference Firn Data
fdm_lat = ncread('./RefData/FGRN11_dens_Apr2012.nc', 'lat');
fdm_lon = ncread('./RefData/FGRN11_dens_Apr2012.nc', 'lon');
fdm_depth = ncread('./RefData/FGRN11_dens_Apr2012.nc', 'depth');
fdm_density = ncread('./RefData/FGRN11_dens_Apr2012.nc', 'dens');

% Load the background firn density and radar calibration constants
load('./RefData/BackgroundDensity.mat');
load('./RefData/2017CalibrationConstants.mat');

% Set the parameters space sweep axes
thickness = 0.01:0.005:0.3;     % Layer thicknesses to test
ice_density = 0.6:0.01:0.92;    % Layer densities to test

% Speed of light in a vacuum (m/s)
c = 299792458;

% Load the radiometric references calculate for a range of background
% densities at the above parameter spacing
rho = 350:1:550;
load('./RefData/RadiometricReferenceCube.mat');

% Flight transects to process
flights = ["20170422_01"];

% Location of layer detection data
data_dir = './DemoData/';
       
for w = 1:length(flights)
    
    fprintf('Processing %s\n', flights(w));
    
    % Load data file and the associated time axis
    data_file = strcat(data_dir, 'Demo_Layer_Detections.mat');
    time_file = strcat(data_dir, 'Demo_TimeRef2.mat');
    load(data_file);
    load(time_file);
    
    extent = length(results);
    
    density_result = cell(1,length(results));
    for k = 1:extent
        
        if k < extent
            fprintf('%d...', k);
        else
            fprintf('%d\n', k);
        end
        
        % Add calibration constant
        data = results{k}.data + calib_constant(10);  % hard coded for demo data
        
        % Set up data structures to save results
        tmp.inversion_max = -1000*ones(size(data));
        tmp.inversion_min = 1000*ones(size(data));
        tmp.lat = results{k}.lat;
        tmp.lon = results{k}.lon;
        tmp.elev = results{k}.elev;
        tmp.base_density = zeros(size(data));
        
        for m = 1:size(data,2)
            
            % Find the nearest background density data points
            opt = abs(fdm_lat - results{k}.lat(m)).^2 + abs(fdm_lon - results{k}.lon(m)).^2;
            [value, index] = min(opt(:));
            [row, col] = ind2sub(size(opt), index);
            
            % Calculate a reference time axis based on the FDM density
            % profile at this location
            n_depth = 0:0.01:20;
            n = 1 + 0.845*(squeeze(fdm_density(row, col, :))/1000);  % Kovacs density to permittivity conversion
            n = vertcat(n(1), n);
            d = vertcat(0, fdm_depth);
            n_interp = interp1(d, n, n_depth);
            true_time = cumsum((2.*diff(n_depth).*n_interp(1:end-1))/c);
            true_time = [0 true_time];
            
            for q = 1:size(data,1)
                if data(q,m) ~= data(1,1)
                    % Convert radargram time axis to depth
                    [~, ind] = min(abs(true_time - time(q)));
                    depth = n_depth(ind);
                    
                    % Find the background dry firn density at this location
                    % and depth
                    [~, ind] = min(abs(background_depth - depth));
                    bg_density = background_density(row, col, ind);
                    
                    tmp.base_density(q,m) = bg_density;
                    
                    % Find the appropriate radiometric reference for that
                    % background density
                    [~, bg_ind] = min(abs(rho - bg_density));
                    radref = 10*log10(squeeze(radrefcube(:,:,bg_ind)));
                    
                    upper = max(max(radref));
                    lower = min(min(radref));
                    
                    % If the observations are within the range predicted by the
                    % radiometric reference
                    if data(q,m) <= upper && data(q,m) >= lower
                        % Find the contour of possible solutions in layer
                        % thickness/layer density space
                        M = contourc(ice_density, thickness, radref, [data(q,m) data(q,m)]);
                        [x,y] = C2xyz(M);
                        % Find the minimum and maximum densities
                        % across all those possible solutions for each
                        % solution contour
                        for p = 1:length(x)
                            if min(x{p}) < tmp.inversion_min(q,m) && length(x{p}) > 10
                                tmp.inversion_min(q,m) = min(x{p});
                            end
                            if max(x{p}) > tmp.inversion_max(q,m) && length(x{p}) > 10
                                tmp.inversion_max(q,m) = max(x{p});
                            end
                        end
                        % Record the absolute highest and lowest possible
                        % solutions for the density
                        if tmp.inversion_min(q,m) == 1000
                            for p = 1:length(x)
                                if min(x{p}) < tmp.inversion_min(q,m)
                                    tmp.inversion_min(q,m) = min(x{p});
                                end
                                if max(x{p}) > tmp.inversion_max(q,m)
                                    tmp.inversion_max(q,m) = max(x{p});
                                end
                            end
                        end
                        if tmp.inversion_max(q,m) == -1000
                            for p = 1:length(x)
                                if min(x{p}) < tmp.inversion_min(q,m)
                                    tmp.inversion_min(q,m) = min(x{p});
                                end
                                if max(x{p}) > tmp.inversion_max(q,m)
                                    tmp.inversion_max(q,m) = max(x{p});
                                end
                            end
                        end
                    % If observations exceed the power predicted by
                    % radiometric reference, assume the layer is definitely
                    % ice
                    elseif data(q,m) > upper
                        tmp.inversion_min(q,m) = 0.92;
                        tmp.inversion_max(q,m) = 0.92;
                    % If observations are less than the predicted
                    % radiometric reference, assume that layer is not ice
                    % and assign it as a non-detection
                    elseif data(q,m) < lower
                        tmp.inversion_min(q,m) = NaN;
                        tmp.inversion_max(q,m) = NaN;
                    end
                else
                    % If no detection was made, then inverted density is
                    % NaN
                    tmp.inversion_min(q,m) = NaN;
                    tmp.inversion_max(q,m) = NaN;
                end
            end
        end
        
        density_result{k} = tmp;
        % Save intermediate files for debugging
        if mod(k,5) == 0
            out_file = strcat(data_dir, flights(w), '_DensityTemp.mat');
            save(out_file, 'density_result');
        end
        
    end
    
    % Save final results
    out_file = strcat('./DerivedData/Demo_Density.mat');
    save(out_file, 'density_result');
    
end
