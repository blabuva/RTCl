function [paramvals] = change_params(namesToChange, valsToChange, paramnames, paramvals, varargin)
%% Change parameter values
% Usage: [paramvals] = change_params(namesToChange, valsToChange, paramnames, paramvals, varargin)
% Outputs:
%       paramvals   - updated paramvals
% Arguments:     
%       namesToChange   - the name(s) of the parameter(s) to change
%                   must be a string/char vec or 
%                       a cell array of strings/char vecs
%       valsToChange    - the values(s) of the parameter(s) to change
%                   must be a numeric array
%       paramnames  - a cell array of all parameter names
%                   must be a cell array of strings or character arrays
%       paramvals   - a numeric array of all parameter values
%                   must be a numeric array of same length as paramnames
%       varargin    - 'ExperimentName': name of the experiment of interest
%                   must be an unambiguous, case-insensitive match 
%                       to one of the following: 
%                       'RTCl'  - chloride accumulation/extrusion model for 
%                                 Peter; updates 'cp_amp', 'cp_per', 'cp_num'
%                       'm3ha'  - GABA-B receptor network
%                       'noexp' - no experiment provided; do nothing
%                   default == 'noexp'
%
% Requires:    
%        /home/Matlab/Adams_Functions/find_ind_str_in_cell.m
%        /home/Matlab/Adams_Functions/update_params.m
%
% Used by:    
%        /media/adamX/RTCl/neuronlaunch.m
%        /media/adamX/m3ha/m3ha_launch.m
%
% 2017-05-03 Moved from update_params.m
% 2017-05-03 Moved trial number update back to neuronlaunch.m
% 2017-05-03 Changed so that it reads namesToChange and valsToChange
% 2018-05-08 Changed tabs to spaces and limited width to 80

%% Hard-coded parameters
validExperiments = {'RTCl', 'm3ha', 'noexp'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 4
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to an input Parser
addRequired(iP, 'namesToChange', ...
    @(x) validateattributes(x, {'char', 'string', 'cell'}, {'nonempty'}));
addRequired(iP, 'valsToChange', ...
    @(x) validateattributes(x, {'numeric'}, {'nonempty'}));
addRequired(iP, 'paramnames', ...
    @(x) assert(iscell(x) && ...
                (min(cellfun(@ischar, x)) || min(cellfun(@isstring, x))), ...
        'Second input must be a cell array of strings or character arrays!'));
addRequired(iP, 'paramvals', ...
    @(x) validateattributes(x, {'numeric', 'logical'}, {'vector'}));

% Add parameter-value pairs to the input Parser
addParameter(iP, 'ExperimentName', 'noexp', ...
    @(x) any(validatestring(x, validExperiments)));

% Read from the input Parser
parse(iP, namesToChange, valsToChange, paramnames, paramvals, varargin{:});
experimentname = validatestring(iP.Results.ExperimentName, validExperiments);

% Check relationships between arguments
if length(paramnames) ~= length(paramvals)
    error('length of paramnames and paramvals are unequal!!');
end

%% Change parameter(s)
if iscell(namesToChange)
    % One or more parameters to change
    ntochange = numel(namesToChange);    % number of parameters to change
    for p = 1:ntochange
        indp = find_ind_str_in_cell(namesToChange{p}, paramnames);
        paramvals(indp) = valsToChange(p);
    end
else
    % Only one parameter to change
    indp = find_ind_str_in_cell(namesToChange, paramnames);
    paramvals(indp) = valsToChange;
end

%% Update dependent parameters for particular experiments
[paramvals] = update_params(paramnames, paramvals, ...
                            'ExperimentName', experimentname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:                                          

%}
