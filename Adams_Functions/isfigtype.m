function [results, figtypes] = isfigtype (strings, varargin)
%% Check whether a string or each string in a cell array is a valid figure type
% Usage: [results, figtypes] = isfigtype (strings, varargin)
% Outputs:    
%       results     - indication of whether the specified name is a figure type
%                   specified as a logical array
%       figtypes    - validated figtypes, if any
%                   specified as a string/char-vec or 
%                       a cell array of strings/char-vecs
%                   returns the shortest match if matchmode == 'substring' 
%                       (sames as validatestring())
% Arguments:        
%       strings     - string or strings to check
%                   must be a string/char-vec or 
%                       a cell array of strings/char-vecs
%       varargin    - 'ValidateMode': whether to validate string and 
%                       throw error if string is not a substring of a figtype
%                   must be logical 1 (true) or 0 (false)
%                   default == false
%                   - 'MatchMode': the matching mode
%                   must be an unambiguous, case-insensitive match 
%                       to one of the following:
%                       'exact'         - string must be exact
%                       'substring'     - string can be a substring
%                   if 'ValidateMode' is 'true', matching mode is 
%                       automatically 'substring'
%                   default == 'substring'
%
% Requires:
%       /home/Matlab/Adams_Functions/validate_string.m
%
% Used by:    
%       /home/Matlab/function_template.m
%       /home/Matlab/minEASE/detect_gapfree_events.m
%       /home/Matlab/Adams_Functions/plot_tuning_curve.m
%       /home/Matlab/Adams_Functions/plot_tuning_map.m
%       /media/adamX/RTCl/raster_plot.m
%       /media/adamX/RTCl/single_neuron.m
%       /media/adamX/RTCl/tuning_curves.m
%       /media/adamX/RTCl/tuning_maps.m
%
% File History:
% 2017-05-09 Created by Adam Lu
% 2017-05-23 Modified linewidth and indentation
% 2018-05-08 Changed tabs to spaces and limited width to 80
% 

%% Taken from Matlab 2017a Documentation
possible_figtypes = {'png', 'fig', 'm', 'mfig', ...
        'jpeg', 'tiff', 'tiffn', 'meta', ...
        'bmpmono', 'bmp', 'bmp16m', 'bmp256', ...
        'hdf', 'pbm', 'pbmraw', 'pcxmono', ...
        'pcx24b', 'pcx256', 'pcx16', 'pgm', ...
        'pgmraw', 'ppm', 'ppmraw'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to an Input Parser
addRequired(iP, 'strings', ...                  % string or strings to check
    @(x) assert(ischar(x) || ...
                iscell(x) && (min(cellfun(@ischar, x)) || ...
                min(cellfun(@isstring, x))) || isstring(x) , ...
                ['strings must be either a string/character array ', ...
                'or a cell array of strings/character arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ValidateMode', false, ...     % whether to validate string
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'MatchMode', 'substring', ...  % the matching mode
    @(x) any(validatestring(x, {'exact', 'substring'})));

% Read from the Input Parser
parse(iP, strings, varargin{:});
validatemode = iP.Results.ValidateMode;
matchmode = iP.Results.MatchMode;

% Display warning message if some inputs are unmatched
if ~isempty(fieldnames(iP.Unmatched))
    fprintf('WARNING: The following name-value pairs could not be parsed: \n');
    disp(iP.Unmatched);
end

%% Check strings and validate
if iscell(strings)
    figtypes = cellfun(@(x) validate_string(x, possible_figtypes, ...
                                            'ValidateMode', validatemode, ...
                                            'MatchMode', matchmode), ...
                                            strings, 'UniformOutput', false);
    results = ~cellfun(@isempty, figtypes);
elseif ischar(strings)
    figtypes = validate_string(strings, possible_figtypes, ...
                                'ValidateMode', validatemode, ...
                                'MatchMode', matchmode);
    results = ~isempty(figtypes);
else
    error(['input argument is in the wrong format! ', ...
            'Type ''help isfigtype'' for usage']);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}
