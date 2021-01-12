% Author: Riley Culberg
% Date: 12/7/2020
%
% This function finds the crossover points on a series of flight transects.
% This function does not return self-matches - only crossovers between
% different transects.
%
% Inputs:
% survey_nums - vector of the survey IDs, must be the same length as
% flight_tracks
% flight_tracks - cell array of Polar Stereographic coordinates for each
% radar trace
%   flight_tracks{m}.x - vector of x coordinates
%   flights_track{m}.y - vector of y coordinates
% tolerance - scalar that sets the maximum distance in meters between two
% points such at it is still considered to be a crossover
% num_xover_guess - scalar, estimated number of crossovers to allow
% preallocation of the output cell array, should be set to some value
% larger than the true number of crossovers
% mean_pow_radius - scalar, radius in meters around a point inside which
% other crossovers should be excluded as duplicates
%
% Ouputs:
% xovers - cell array of crossover points
%   xovers{m}.survey1 - ID of first survey in crossover
%   xovers{m}.survey2 - ID of first survey in crossover
%   xovers{m}.survey1_x - x coordinate of crossover point in first survey
%   xovers{m}.survey1_y - y coordinate of crossover point in first survey
%   xovers{m}.survey2_x - x coordinate of crossover point in second survey
%   xovers{m}.survey2_y - y coordinate of crossover point in second survey
% -------------------------------------------------------------------------

function xovers = find_xover(survey_nums, flight_tracks, tolerance, num_xover_guess, mean_pow_radius)

    pairings = nchoosek(survey_nums, 2);
    xovers = cell(1, num_xover_guess);
    
    count = 0;
    for k = 1:size(pairings,1)
        fprintf('Calculating Survey %d Crossovers with Survey %d\n', pairings(k,1), pairings(k,2));
        for m = 1:length(flight_tracks{pairings(k,1)}.x)
            [min_dist, ind] = min(abs(flight_tracks{pairings(k,2)}.x - flight_tracks{pairings(k,1)}.x(m)).^2 + ...
                           abs(flight_tracks{pairings(k,2)}.y - flight_tracks{pairings(k,1)}.y(m)).^2 );
            if sqrt(min_dist) <= tolerance
                count = count + 1;
                location.survey1 = pairings(k,1);
                location.survey2 = pairings(k,2);
                location.survey1_x = flight_tracks{pairings(k,1)}.x(m);
                location.survey1_y = flight_tracks{pairings(k,1)}.y(m);
                location.survey2_x = flight_tracks{pairings(k,2)}.x(ind);
                location.survey2_y = flight_tracks{pairings(k,2)}.y(ind);
                xovers{count} = location;
            end
        end
    end
    
    % Delete extra cells
    full_cells = find(cellfun(@isempty,xovers));
    if ~isempty(full_cells)
        xovers(full_cells) = [];
    end
    
    % Remove duplicated crossovers that occur within the proscribed radius
    % by finding all adjacent crossovers and only keeping the one with the
    % minimum distance between the survey1 and survey2 points
    count = 1;
    while count < length(xovers)
        index = [];
        for m = count:length(xovers)
            dist = sqrt(abs(xovers{count}.survey1_x - xovers{m}.survey1_x).^2 + abs(xovers{count}.survey1_y - xovers{m}.survey1_y).^2);
            if dist <= mean_pow_radius
                index = [index m];
            end
        end
        distance = zeros(size(index));
        
        for n = 1:length(index)
            distance(n) = sqrt(abs(xovers{index(n)}.survey1_x - xovers{index(n)}.survey2_x).^2 + abs(xovers{index(n)}.survey1_y - xovers{index(n)}.survey2_y).^2);
        end
        [~, ind] = min(distance);
        index(ind) = [];
        xovers(index) = [];
        count = count + 1;
    end

end