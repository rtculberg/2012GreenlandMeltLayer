% Author: Riley Culberg
% Date: 12/7/2020
%
% This function finds the exact crossover points for a series of flight
% transects given approximate crossover locations to begin the search. 
%
% Inputs:
% flight_tracks - cell array of Polar Stereographic coordinates for each
% radar trace
%   flight_tracks{m}.x - vector of x coordinates
%   flights_track{m}.y - vector of y coordinates
% xovers - cell array of crossover points
%   xovers{m}.survey1 - ID of first survey in crossover
%   xovers{m}.survey2 - ID of first survey in crossover
%   xovers{m}.survey1_x - x coordinate of crossover point in first survey
%   xovers{m}.survey1_y - y coordinate of crossover point in first survey
%   xovers{m}.survey2_x - x coordinate of crossover point in second survey
%   xovers{m}.survey2_y - y coordinate of crossover point in second survey
% search area - radius around the approximate crossovers to search in
% number of traces (should be the same amount as was used to downsample
% flight tracks for initial rough pass or slightly more) (scalar double)
% rf - radius over which to average the surface power around the crossover
% point in number of traces (scalar double)
% 
% Ouputs:
% fine_xovers - cell array of crossover points
%   fine_xovers{m}.survey1 - ID of first survey in crossover
%   fine_xovers{m}.survey2 - ID of first survey in crossover
%   fine_xovers{m}.survey1_x - x coordinate of crossover point in first survey
%   fine_xovers{m}.survey1_y - y coordinate of crossover point in first survey
%   fine_xovers{m}.survey2_x - x coordinate of crossover point in second survey
%   fine_xovers{m}.survey2_y - y coordinate of crossover point in second survey 
% -------------------------------------------------------------------------

function fine_xovers = find_fine_xover(flight_tracks, xovers, search_area, rf)

    fine_xovers = cell(1, length(xovers));
    
    for k = 1:length(xovers)
        
        % Apprximate x,y location of crossover
        s1 = xovers{k}.survey1;
        s2 = xovers{k}.survey2;
        
        % Find the trace indices of the approximate crossover point
        [~, ind1] = min(abs(flight_tracks{s1}.x - xovers{k}.survey1_x).^2 + abs(flight_tracks{s1}.y - xovers{k}.survey1_y).^2);
        [~, ind2] = min(abs(flight_tracks{s2}.x - xovers{k}.survey2_x).^2 + abs(flight_tracks{s2}.y - xovers{k}.survey2_y).^2);
        
        % Define the search area, making sure we can't go past the start of
        % end of a flight track
        if ind1+search_area > length(flight_tracks{s1}.x)
            search_area = length(flight_tracks{s1}.x) - ind1;
        elseif ind1-search_area < 0
            search_area = ind1 - 1;
        elseif ind1+search_area > length(flight_tracks{s1}.y)
            search_area = length(flight_tracks{s1}.y) - ind1;
        end
        
        if ind2+search_area > length(flight_tracks{s2}.x)
            search_area = length(flight_tracks{s2}.x) - ind2;
        elseif ind2-search_area < 0
            search_area = ind2 - 1;
        elseif ind2+search_area > length(flight_tracks{s2}.y)
            search_area = length(flight_tracks{s2}.y) - ind1;
        end
        
        % Clip flight tracks to the search area
        s1_search_x = flight_tracks{s1}.x(ind1-search_area:ind1+search_area);
        s1_search_y = flight_tracks{s1}.y(ind1-search_area:ind1+search_area);
        s2_search_x = flight_tracks{s2}.x(ind2-search_area:ind2+search_area);
        s2_search_y = flight_tracks{s2}.y(ind2-search_area:ind2+search_area);
        
        % Search for the two points with the minimum distance between the
        % different transects
        min_dist = 30000;
        min_ind_s1 = 0;
        min_ind_s2 = 0;
        for m = 1:length(s1_search_x)
            [M, ind] = min(abs(s2_search_x - s1_search_x(m)).^2 + abs(s2_search_y - s1_search_y(m)).^2);
            if sqrt(M) < min_dist
                min_dist = sqrt(M);
                min_ind_s1 = m;
                min_ind_s2 = ind;
            end
        end
        
        % Recalibrate crossover index for full fligh track
        abs_s1_ind = min_ind_s1 + (ind1 - search_area) - 1;
        abs_s2_ind = min_ind_s2 + (ind2 - search_area) - 1;
        
        % Save location of crossover
        location.survey1 = s1;
        location.survey2 = s2;
        location.survey1_x = s1_search_x(min_ind_s1);
        location.survey1_y = s1_search_y(min_ind_s1);
        location.survey2_x = s2_search_x(min_ind_s2);
        location.survey2_y = s2_search_y(min_ind_s2);
        
        % Setup of averaging radius for surface power making sure that we
        % don't got past the start or end of a flight transect
        if abs_s1_ind + rf > length(flight_tracks{s1}.x)
            rf = length(flight_tracks{s1}.x) - abs_s1_ind;
        elseif abs_s1_ind - rf <= 0
            rf = ind1 - 1;
        elseif abs_s1_ind + rf > length(flight_tracks{s1}.y)
            rf = length(flight_tracks{s1}.y) - abs_s1_ind;
        end
        
        if abs_s2_ind + rf > length(flight_tracks{s2}.x)
            rf = length(flight_tracks{s2}.x) - abs_s2_ind;
            disp(rf);
        elseif abs_s2_ind - rf <= 0
            rf = ind2 - 1;
        elseif abs_s2_ind + rf > length(flight_tracks{s2}.y)
            rf = length(flight_tracks{s2}.y) - abs_s2_ind;
        end
        
        % Average surface power at the crossover
        location.survey1_pow = mean(flight_tracks{s1}.surf_geo_pow(abs_s1_ind - rf:abs_s1_ind + rf));
        location.survey2_pow = mean(flight_tracks{s2}.surf_geo_pow(abs_s2_ind - rf: abs_s2_ind + rf));
        
        fine_xovers{k} = location;
    end

end