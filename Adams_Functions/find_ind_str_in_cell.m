function [indices, elements] = find_ind_str_in_cell(str, cellArray, varargin)
%% Find all indices of a particular string in a cell array
% Usage: [indices, elements] = find_ind_str_in_cell(str, cellArray, varargin)
% Explanation:
%   This works like the strcmp() or strcmpi function in Matlab, 
%       especially when 'SearchMode' == 'exact'.
%   However, this returns indices and not logical arrays,
%       and returns the corresponding elements.
%   Also, default is 'SearchMode' == 'substrings', which allows str to be 
%       a substring of a match in cellArray.
% Outputs:
%       indices     - indices of the cell array containing that exact string
%                       or containing a substring or all substrings provided; 
%                       could be empty
%                   specified as a numeric array
%       elements    - elements of the cell array corresponding to those indices
%                   specified as a cell array if more than one indices 
%                       or the element if only one index; or an empty string
% Arguments:
%       str         - a string or a substring or 
%                       a cell array of substrings of interest
%                   must be a string/character array or 
%                       a cell array of strings/character arrays
%       cellArray   - a cell array that contains strings
%                   must be a cell array of strings/character arrays
%       varargin    - 'SearchMode': the search mode
%                   must be an unambiguous, case-insensitive match to one of:
%                       'exact'         - str must be identical to 
%                                           an element in cellArray
%                       'substrings'    - str can be a substring or 
%                                           a cell array of substrings
%                   if searchMode is 'exact', str cannot be a cell array
%                   default == 'substrings'
%                   - 'IgnoreCase': whether to ignore differences in letter case
%                   must be logical 1 (true) or 0 (false)
%                   default == false
%                   - 'MaxNum': maximum number of indices to find
%                   must be a positive integer
%                   default == numel(cellArray)
%
% Requires:
%       /home/Matlab/Adams_Functions/intersectm.m
%
% Used by:    
%       /media/adamX/m3ha/data_dclamp/dclampPassiveFitter.m
%       /media/adamX/m3ha/data_dclamp/PlotHistogramsRefineThreshold.m
%       /media/adamX/m3ha/data_dclamp/test_sweep.m
%       /media/adamX/m3ha/data_dclamp/remove_E092810_0000.m
%       /media/adamX/m3ha/data_dclamp/compare_statistics.m
%       /media/adamX/m3ha/optimizer4gabab/optimizergui_4compgabab.m
%       /media/adamX/m3ha/optimizer4gabab/optimizer_4compgabab.m
%       /media/adamX/m3ha/optimizer4gabab/compare_neuronparams.m
%       /media/adamX/RTCl/neuronlaunch.m
%       /media/adamX/RTCl/raster_plot.m
%       /media/adamX/RTCl/tuning_curves.m
%       /media/adamX/RTCl/single_neuron.m
%       /home/Matlab/Adams_Functions/validate_string.m
%       /home/Matlab/Adams_Functions/update_params.m
%       /home/Matlab/Adams_Functions/increment_editbox.m
%       /home/Matlab/minEASE/read_params.m
%       /home/Matlab/minEASE/gui_examine_events.m
%       /home/Matlab/EEG_gui/EEG_gui.m
%
% 2016-09--- Created
% 2016-10-13 moved to Adams_Functions
% 2016-11-30 Added searchMode
% 2017-04-05 Fixed the size of str_cell so that it can take column or row arrays
% 2017-04-26 Now str can be a cell array of substrings too
% 2017-04-27 Improved inputParser scheme
% 2017-05-09 Added elements as output
% 2017-05-25 Changed line width and indentation
% 2017-06-09 Fixed the returned element to be of original case
% 2018-05-01 Added MaxNum as a parameter

%% Default values for optional arguments
searchModeDefault = 'substrings';       % default search mode
ignoreCaseDefault = false;              % whether to ignore case by default
maxNumDefault = [];                     % will be changed to numel(cellArray)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'str', ...              % a string/substrings of interest
    @(x) assert(ischar(x) || iscell(x) && (min(cellfun(@ischar, x)) || ...
                min(cellfun(@isstring, x))) || isstring(x), ...
                ['First input must be either a string/character array ', ...
                'or a cell array of strings/character arrays!']));
addRequired(iP, 'cellArray', ...        % a cell array that contains strings
    @(x) assert(iscell(x) && (min(cellfun(@ischar, x)) || ...
                min(cellfun(@isstring, x))), ...
                ['Second input must be a cell array ', ...
                'of strings/character arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SearchMode', searchModeDefault, ...   % the search mode
    @(x) any(validatestring(x, {'exact', 'substrings'})));
addParameter(iP, 'IgnoreCase', ignoreCaseDefault, ...   % whether to ignore case
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'MaxNum', maxNumDefault, ...       % maximum number of indices
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer'}));

% Read from the Input Parser
parse(iP, str, cellArray, varargin{:});
searchMode = validatestring(iP.Results.SearchMode, {'exact', 'substrings'});
ignoreCase = iP.Results.IgnoreCase;
maxNum = iP.Results.MaxNum;

% Check relationships between arguments
if strcmpi(searchMode, 'exact') && iscell(str)
    error('First input cannot be a cell array if searchMode is 1!');
end

%% Prepare for the search
% Count the number of indices
nIndices = numel(cellArray);

% Set the maximum number of indices if not provided
if isempty(maxNum)
    maxNum = nIndices;
end

%% Find the indices
if strcmpi(searchMode, 'exact')             % String must be exact

    % Construct a cell array with the same size as cellArray 
    %   but with str duplicated throughout
    str_cell = cell(size(cellArray));    % a cell array to store copies of str
    for k = 1:numel(cellArray)
        % Make the kth element the same as str
        str_cell{k} = str;
    end

    % Find indices that corresponds to str exactly in cellArray, 
    %   case-insensitive if IgnoreCase is set to true
    if ignoreCase
        indices = find(cellfun(@strcmpi, cellArray, str_cell), maxNum);
    else
        indices = find(cellfun(@strcmp, cellArray, str_cell), maxNum);
    end

elseif strcmpi(searchMode, 'substrings')    % String can be a substring 
                                            % or a cell array of substrings

    % Convert each string to lower case if IgnoreCase is set to true
    if ignoreCase
        cellArrayMod = cellfun(@lower, cellArray, 'UniformOutput', false);
    else
        cellArrayMod = cellArray;
    end

    % Find indices that contain str in cellArrayMod
    if iscell(str)        % if str is a cell array of substrings
        % Convert each substring to lower case if IgnoreCase is set to true
        if ignoreCase
            strMod = cellfun(@lower, str, 'UniformOutput', false);
        else
            strMod = str;
        end

        % Find the indices that contain each substring
        numstrs = numel(strMod);
        indices_eachstr = cell(1, numstrs);
        for k = 1:numstrs
            indicesarray = strfind(cellArrayMod, strMod(k));    
            indices_eachstr{k} = find(~cellfun(@isempty, indicesarray));
        end

        % Find the indices that contain all substrings by intersection
        indices = intersectm(indices_eachstr);

        % If more than maxNum indices found, 
        %   restrict to the first maxNum indices
        if length(indices) > maxNum
            indices = indices(1:maxNum);
        end
    else                    % if str is a single substring
        % Convert substring to lower case if IgnoreCase is set to true
        if ignoreCase
            strMod = lower(str);
        else
            strMod = str;
        end

        % Find the indices that contain the substring
        indicesarray = strfind(cellArrayMod, strMod);
        indices = find(~cellfun(@isempty, indicesarray), maxNum);
    end

end

%% Return the elements too
if ~isempty(indices) && numel(indices) > 1
    elements = cellArray(indices);
elseif ~isempty(indices)
    elements = cellArray{indices};
else
    elements = '';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}
