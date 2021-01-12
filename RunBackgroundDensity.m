% Author: Riley Culberg
% Date: 12/7/2020
%
% This script calculates the mean background firn density as a function of 
% depth from the IMAU FDM firn model outputs for all of Greenland.
%
% Outputs:
%   background_depth - depth axis in meters (Mx1 double vector)
%   background_density - background density in g/cm^3 for all of Greenland
%   (NxPxM double vector). The N axis is latitude and P is longitude as
%   defined by the IMAU FDM grid.
%
% Data dependencies:
%   - IMAU FDM data - ./RefData/FGRN11_dens_Apr2012.nc
% ------------------------------------------------------------------------

clear;

addpath(genpath('./ReferenceFunctions/'));

% Load FDM firn model data
fdm_lat = ncread('./RefData/FGRN11_dens_Apr2012.nc', 'lat');
fdm_lon = ncread('./RefData/FGRN11_dens_Apr2012.nc', 'lon');
fdm_depth = ncread('./RefData/FGRN11_dens_Apr2012.nc', 'depth');
fdm_density = ncread('./RefData/FGRN11_dens_Apr2012.nc', 'dens');

% Set up the background depth axis
background_depth = 0:0.1:20;
background_density = zeros(size(fdm_density,1), size(fdm_density,2), length(background_depth));
for k = 1:size(fdm_density,1)
    for m = 1:size(fdm_density,2)
        
        line = squeeze(fdm_density(k,m,:));
        d = fdm_depth;
        
        if ~isnan(line(1))
            
            % Find the lower envelope of the density profile so
            % that we are not fitting ice layers, but only the background
            % dry firn density profile
            
            x = gradient(gradient(line));
            ind = find(x < 0);
            line(ind) = [];
            d(ind) = [];
            
            x = gradient(gradient(line));
            ind = find(x < 0);
            line(ind) = [];
            d(ind) = [];
            
            ind = find(~isnan(line));
            line = line(ind);
            d = d(ind);
            
            cutoff = 60;
            [~, stop] = min(abs(d - cutoff));
            if stop > length(line)
                stop = length(line);
            end
            
            % Fit a biexponential function to the lower envelope
            y = fit(d(1:stop), line(1:stop), 'exp2');
            model = y.a*exp(y.b*background_depth) + y.c*exp(y.d*background_depth);
            
            background_density(k,m,:) = model;
            
        else
            background_density(k,m,:) = NaN;
        end
        
    end
end

% Save results
save('./DerivedData/BackgroundDensity.mat', 'background_depth', 'background_density');
