function [pchnames, pchvalues, ntrials, nump, pvalues, nperp] = make_loopedparams (loopmode, pnames, plabels, pislog, pmin, pmax, pinc, varargin)
%% Construct parameters to change for each trial from loopmode, pnames, pislog, pmin, pmax, pinc 
% Usage: [pchnames, pchvalues, ntrials, nump, pvalues, nperp] = make_loopedparams (loopmode, pnames, plabels, pislog, pmin, pmax, pinc, varargin)
% Outputs:    
%       pchnames    - a cell array of parameter names or ordered pairs 
%                       of parameter names for each trial
%       pchvalues   - a numeric array of parameter values or a cell array 
%                       of ordered paris of parameter values for each trial
%       ntrials     - total number of trials
%       nump        - number of different parameters
%       pvalues     - a cell array of arrays of parameter values
%       nperp       - number of parameter values for each parameter
% Arguments:    
%       loopmode    - how to loop through parameters: 'cross' or 'grid'
%                   must be an unambiguous, case-insensitive match 
%                       to one of the following: 
%                       'cross' - Loop through each parameter 
%                                   while fixing others
%                       'grid'  - Loop through all possible 
%                                   combinations of parameters
%       pnames      - names of parameters to loop through
%                   must be a cell array of strings or character arrays
%       plabels     - labels of parameters to loop through
%                   must be a cell array of strings or character arrays
%       pislog      - whether increments of parameters is in log
%                   must be a vector of logicals or 0s/1s
%       pmin        - minimum values of parameters to loop through
%                   must be a numeric vector
%       pmax        - maximum values of parameters to loop through
%                   must be a numeric vector
%       pinc        - increments of parameters to loop through
%                   must be a numeric vector
%       varargin    - 'OutFolder': directory to place outputs, 
%                                   e.g. '/media/shareX/share/'
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'FileLabel': label for output files
%                   must be a string scalar or a character vector
%                   default == datestr(clock, 30);
%                   - 'NCells': number of cells
%                   must be a positive integer
%                   default == 0;
%                   - 'ActMode': activation mode
%                   must be a nonegative integer
%                   default == 0;
%
% Requires:
%       cd/all_ordered_pairs.m
%
% Used by:    
%       /media/adamX/RTCl/neuronlaunch.m
%
% File History:
% 2017-05-03 Moved from /media/adamX/RTCl/neuronlaunch.m
% 2018-05-08 Changed tabs to spaces and limited width to 80
% 

%% Hard-coded parameters
possibleLoopmodes = {'cross', 'grid'};
loopedparamsfile = 'loopedparams.mat';
            % must be consistent with get_loopparams.m & combine_loopparams.m

%% Default values for optional arguments
defaultFileLabel = datestr(clock, 30);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 7
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to an input Parser
addRequired(iP, 'loopmode', ... % how to loop through parameters
    @(x) any(validatestring(x, possibleLoopmodes)));
addRequired(iP, 'pnames', ...   % names of parameters to loop through
    @(x) assert(iscell(x) && ...
                (min(cellfun(@ischar, x)) || min(cellfun(@isstring, x))), ...
        'Second input must be a cell array of strings or character arrays!'));
addRequired(iP, 'plabels', ...  % labels of parameters to loop through
    @(x) assert(iscell(x) && ...
                (min(cellfun(@ischar, x)) || min(cellfun(@isstring, x))), ...
        'Second input must be a cell array of strings or character arrays!'));
addRequired(iP, 'pislog', ...   % whether increments of parameters is in log
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary', 'vector'}));
addRequired(iP, 'pmin', ...     % minimum values of parameters to loop through
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'pmax', ...     % maximum values of parameters to loop through
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'pinc', ...     % increments of parameters to loop through
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Add parameter-value pairs to the input Parser
addParameter(iP, 'OutFolder', '', ...   % directory to place outputs
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
%    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));    
                                                    % introduced after R2016B
addParameter(iP, 'FileLabel', defaultFileLabel, ...     % label for output files
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
%    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                    % introduced after R2016B
addParameter(iP, 'NCells', 0, ...       % number of cells
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer'}));
addParameter(iP, 'ActMode', 0, ...      % activation mode
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'integer'}));

% Read from the input Parser
parse(iP, loopmode, pnames, plabels, pislog, pmin, pmax, pinc, varargin{:});
loopmode = validatestring(loopmode, possibleLoopmodes);
outfolder = iP.Results.OutFolder;
filelabel = iP.Results.FileLabel;
ncells = iP.Results.NCells;
actmode = iP.Results.ActMode;

% Set default arguments
if isempty(outfolder)
    outfolder = pwd;
end

% Check relationships between arguments
if numel(unique([numel(pnames), numel(plabels), numel(pislog), ...
                numel(pmin), numel(pmax), numel(pinc)])) ~= 1
    error('pnames, plabels, pislog, pmin, pmax, pinc not all equal length!!');
end

%% Calculate number of parameters to loop through
nump = numel(pnames);
pvalues = cell(1, nump);        % parameter values that will be used
nperp = zeros(1, nump);            % number of trials per parameter type
for p = 1:nump
    % Generate parameter values that will be used
    if pislog(p)
        pvalues{p} = exp(log(pmin(p)):log(pinc(p)):log(pmax(p)))';
                                            % parameters used as a column vector
    else
        pvalues{p} = (pmin(p):pinc(p):pmax(p))';
                                            % parameters used as a column vector
    end
    nperp(p) = length(pvalues{p});
end

%% Find total number of trials
switch loopmode
case 'cross'
    % Total number of trials in 'cross' mode
    ntrials = sum(nperp);
case 'grid'
    % Total number of trials in 'grid' mode
    ntrials = prod(nperp);
end

%% Create pchnames & pchvalues
switch loopmode
case 'cross'
    pchnames = cell(1, ntrials);    % stores parameter name for each trial
    pchvalues = zeros(1, ntrials);  % stores parameter value for each trial
    ct = 0;                         % counts number of trials used
    for p = 1:nump
        % Store parameter names and values for each trial 
        %   (to keep things in order)
        indices = (ct + 1):(ct + nperp(p));     % indices for this parameter
        pchnames(indices) = repmat(pnames(p), 1, length(indices));
                                                % parameter name for each trial
        pchvalues(indices) = (pvalues{p})';     % parameter value for each trial 

        % Update count of number of files used
        ct = ct + nperp(p);
    end
case 'grid'
    % Repeat each parameter name nperp times
    pnames_rep = cell(1, nump);     % stores each parameter name 
                                    %   repeated nperp times
    for p = 1:nump
        pnames_rep{p} = repmat(pnames(p), nperp(p), 1);
    end

    % Construct all possible ordered pairs of values and corresponding names
    pchvalues = all_ordered_pairs(pvalues); 
                                    % a cell array of all ordered pairs
    pchnames = all_ordered_pairs(pnames_rep);
                                    % a cell array of ordered pairs of 
                                    %   corresponding parameter names
end

%% Save looping variables in a mat file named by the date & time
save(fullfile(outfolder, sprintf('%s_%s', filelabel, loopedparamsfile)), ...
    'loopmode', 'ntrials', 'nump', 'nperp', ...
    'pnames', 'plabels', 'pislog', 'pmin', 'pmax', 'pinc', ...
    'pvalues', 'pchnames', 'pchvalues', ...
    'ncells', 'actmode', '-v7.3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

