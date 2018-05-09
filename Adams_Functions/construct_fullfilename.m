function [fullfilename] = construct_fullfilename (filename, varargin)
%% Constructs full file name based on filename and an optional full directory path with optional suffices and/or Name-Value pairs
% Usage: [fullfilename] = construct_fullfilename (filename, varargin)
% Outputs:
%       fullfilename    - the full path to file constructed
% Arguments:
%       filename    - filename with extension but without directory, 
%                       e.g. 'A100110_0008_18.mat'
%                   must be a string scalar or a character vector
%       varargin    - 'Directory': a full directory path, 
%                       e.g. '/media/shareX/share/'
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'Suffices': suffix(ces) to add to filebase
%                   must be a string/character array or a cell array 
%                       of strings/character arrays
%                   default == 'nosuffices'
%                   - 'NameValuePairs': Name-Value pairs that are changed
%                   must be a 2-element cell array whose first element 
%                       is a string/char array or cell array 
%                       and whose second element is a numeric array
%                   default == {'nopairs', NaN}
%        
%
% Requires:
%       cd/construct_suffix.m
% Used by:
%       /media/adamX/RTCl/neuronlaunch.m
% 
% 2017-03-27 Created
% 2017-05-04 Removed ntrials
% 2017-05-04 Added input Parser scheme
% 2017-05-04 Changed paramname, paramvalue -> 'NameValuePairs' 
%                and made it a parameter 
% 2017-05-04 Made 'suffices' a parameter 
% 2018-05-08 Changed tabs to spaces and limited width to 80

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Add required inputs to an input Parser
iP = inputParser;
addRequired(iP, 'filename', ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
%    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                    % introduced after R2016B

% Add parameter-value pairs to the input Parser
addParameter(iP, 'Directory', '', ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
%    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                    % introduced after R2016B
addParameter(iP, 'Suffices', '', ...
    @(x) assert(ischar(x) || iscell(x) && (min(cellfun(@ischar, x)) || ...
                min(cellfun(@isstring, x))) || isstring(x), ...
                ['Suffices must be either a string/character array ', ...
                    'or a cell array of strings/character arrays!']));
addParameter(iP, 'NameValuePairs', {'', NaN}, ...
    @(x) assert(iscell(x) && numel(x) == 2 ...
            && (ischar(x{1}) || iscell(x{1}) || isstring(x{1})) ...
            && isnumeric(x{2}), ...
        ['NameValuePairs must be a 2-element cell array whose ', ...
            'first element is a string/char array or cell array ', ...
            'and whose second element is a numeric array!']));

% Read from the input Parser
parse(iP, filename, varargin{:});
directory = iP.Results.Directory;
suffices = iP.Results.Suffices;
namevaluepairs = iP.Results.NameValuePairs;

% Set default arguments
if isempty(directory)
    directory = pwd;        % default path is current directory
end

%% Do the job
% Separate filebase and extension from filename
[~, filebase, extension] = fileparts(filename);

% Construct final suffix
finalSuffix = construct_suffix('Suffices', suffices, ...
                                'NameValuePairs', namevaluepairs);

% Construct full file name
if isempty(finalSuffix)                % if nothing provided
    % Construct path based on directory, filebase and extension
    fullfilename = fullfile(directory, ...
                            [filebase, extension]);    
elseif ~isempty(finalSuffix)              % suffix(ces) is(are) provided
    % Construct path based on directory, filebase, final suffix and extension
    fullfilename = fullfile(directory, ...
                            [filebase, '_', finalSuffix, extension]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}
