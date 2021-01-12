% Author: Riley Culberg
% Date: 12/7/2020
%
% This scripts retracks the ice surface in the CReSIS L1B radargrams. 
%
% Outputs:
%   surf - for each flight transect, save a .mat file containing the surf
%   variable which is a 1xM cell array where M is the number of radargram
%   segments in the flight transect. Each cell contains a 1xN vector of
%   integers with the fast time sample location of the ice surface at that
%   particular trace.
% 
% Data dependencies:
% CReSIS L1B radargrams - ./DemoData/Data_20170422_01_158.mat for demoing
% the functionality
% ------------------------------------------------------------------------

clear;

addpath(genpath('./ReferenceFunctions/'));

% List of flight transects and radargram segments to process
flights = {'20170422_01'};
segments{1} = [158];

starter = 1000;               % Index guard to avoid antenna feedthrough
base_dir = './DemoData/';     % Data directory
out_dir = './DerivedData/';   % Location to save output
for k = 1:length(flights)
    fprintf('Date: %s\n', flights{k});
    surf = cell(1, length(segments{k}));
    for m = 1:length(segments{k})
        if m < length(segments{k})
            fprintf('%d...\n', segments{k}(m));
        else
            fprintf('%d\n', segments{k}(m));
        end
        flight_results.surface = [];
        
        % Load Data
        if segments{k}(m) < 10
            file1 = strcat(base_dir, '/Data_', flights{k}, '_00', num2str(segments{k}(m)), '.mat');
            file2 = strcat(base_dir, '/Data_img_01_', flights{k}, '_00', num2str(segments{k}(m)), '.mat');
        elseif segments{k}(m) < 100
            file1 = strcat(base_dir, '/Data_', flights{k}, '_0', num2str(segments{k}(m)), '.mat');
            file2 = strcat(base_dir, '/Data_img_01_', flights{k}, '_0', num2str(segments{k}(m)), '.mat');
        else
            file1 = strcat(base_dir, '/Data_', flights{k}, '_', num2str(segments{k}(m)), '.mat');
            file2 = strcat(base_dir, '/Data_img_01_', flights{k}, '_', num2str(segments{k}(m)), '.mat');
        end
        
        load(file1);

        % Load the initial CRESIS tracked surface
        surf_ind_init = zeros(size(Surface));
        for q = 1:length(Surface)
            [~, surf_ind_init(q)] = min(abs(Time - Surface(q)));
        end
        
        % Setup the tracking window around the surface based on the inital
        % CReSIS picks
        start = NaN*ones(size(surf_ind_init));
        stop = NaN*ones(size(surf_ind_init));
        num = 40;
        for p = 1:length(surf_ind_init)
            if surf_ind_init(p) > starter && surf_ind_init(p) < (size(Data,1) - num)
                start(p) = surf_ind_init(p) - num;
                stop(p) = surf_ind_init(p) + 5;
            elseif surf_ind_init(p) < starter && surf_ind_init(p) < (size(Data,1) - num)
                start(p) = max(surf_ind_init) - num;
                stop(p) = max(surf_ind_init) + 5;
            end
        end
        
        % Find the maximum gradient in the window around the surface
        surf_ind = NaN*ones(1, size(Data,2));
        for p = 1:length(surf_ind)
            if ~isnan(start(k)) && ~isnan(stop(p))
                max_del_y = 0;
                surface = 0;
                for q = start(p):1:stop(p)
                    del_y = 10*log10(Data(q,p)) - 10*log10(Data(q-1,p));
                    if del_y > max_del_y
                        max_del_y = del_y;
                        surface = q;
                    end
                end
                surf_ind(p) = surface;
            end
        end
        
        % Smooth the tracking horizontally to avoid large jumps from
        % surface
        for p = 2:length(surf_ind)-1
            if abs(surf_ind(p) - surf_ind(p-1)) > 3 && abs(surf_ind(p) - surf_ind(p+1)) > 3 && ~isnan(surf_ind(p)) && ~isnan(surf_ind(p-1)) && ~isnan(surf_ind(p+1))
                surf_ind(p) = surf_ind(p-1);
            end
        end
        
        % Find the first maximum after the max gradient and save the index
        % and surface power
        for p = 1:length(surf_ind)
            if ~isnan(surf_ind(p))
                surface = surf_ind(p);
                surf_power = 10*log10(Data(surface,p));
                for q = surface+1:1:(surface + 4)
                    if 10*log10(Data(q,p)) > surf_power
                        surf_power = 10*log10(Data(q,p));
                        surface = q;
                    end
                    if 10*log10(Data(q+1,p)) < 10*log10(Data(q,p))
                        break;
                    end
                end
                surf_ind(p) = surface;
            end
        end
        surf{m} = surf_ind;
    end
    % Save results
    save(strcat(out_dir, 'Demo_Surface.mat'), 'surf'); 
end


