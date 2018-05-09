function intersection = intersectm (arrays)
%% Find the intersection of multiple arrays
% Usage: intersection = intersectm (arrays)
% Outputs:    
%       intersection    - intersection of input arrays
% Arguments:
%       arrays      - a cell array of arrays that will be intersected;
%                       if just an array, return the array
%                   must be a cell array of input arrays that 
%                       can be recognized by the built-in intersect function
% Used by:    
%       /home/Matlab/Adams_Functions/find_ind_str_in_cell.m
%
% 2018-01-10 Created
% 2018-05-08 Fixed the case if arrays is not a cell array
% 2018-05-08 Changed tabs to spaces and limited width to 80


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

if isempty(arrays)
    error('Argument cannot be empty!');
end

%% Return the array if not a cell array
if ~iscell(arrays)
    intersection = arrays;
    return;
end

%% Return the intersection of the contents of the cell array
numarrays = numel(arrays);        % number of arrays
intersection = arrays{1};
if numarrays > 1
    for k = 2:numarrays
        intersection = intersect(intersection, arrays{k});
    end
end