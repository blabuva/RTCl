function [paramvals] = update_params(paramnames, paramvals, varargin)
%% Update dependent parameters for particular experiments
% Usage: [paramvals] = update_params(paramnames, paramvals, varargin)
% Outputs:     paramvals    - updated paramvals
% Arguments:
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
%
% Used by:    
%        /home/Matlab/Adams_Functions/change_params.m
%
% 2017-03-30 Created by Adam Lu
% 2017-05-03 Moved to /media/adamX/RTCl/
% 2017-05-03 Moved parts to change_params.m
% 2017-05-03 Added inputParser scheme
% 2017-05-03 Added 'ExperimentName' as a parameter
% 2017-11-07 Fixed TCgabaaGmax for bicuculline mode
% 2018-05-08 Changed tabs to spaces and limited width to 80
% 2018-05-09 Fixed cpamp calculation for 'RTCl'

%% Hard-coded parameters
validExperiments = {'RTCl', 'm3ha', 'noexp'};

%% m3ha: Template parameters used for GABA-B IPSC conductance waveforms
% Note: amp in Christine's thesis was actually for 200 % G incr
%       This must be consistent with /media/adamX/m3ha/data_dclamp/CountSweeps.m
gtemp_gmax = [4.48; 6.72; 9.76; 4.48]./1000;    % maximal conductance (uS)
gtemp_amp = [16.00; 24.00; 8.88; 6.32]./1000;   % conductance amplitude (uS)
gtemp_Trise = [52.00; 52.00; 38.63; 39.88]; % rising phase time constant (ms)
gtemp_TfallFast = [90.10; 90.10; 273.40; 65.80];% fast decay time constant (ms)
gtemp_TfallSlow = [1073.20; 1073.20; 1022.00; 2600.00];
                                                % slow decay time constant (ms)
gtemp_w = [0.952; 0.952; 0.775; 0.629];         % weight of the fast decay

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

% Add required inputs to an input Parser
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
parse(iP, paramnames, paramvals, varargin{:});
experimentname = validatestring(iP.Results.ExperimentName, validExperiments);

% Check relationships between arguments
if length(paramnames) ~= length(paramvals)
    error('length of paramnames and paramvals are unequal!!');
end

%% Update dependent parameters, part I
switch experimentname
case {'RTCl', 'm3ha'}
    % Find indices of parameters that influence other parameters
    indREdia = find_ind_str_in_cell('REdiam', paramnames);
    indstimf = find_ind_str_in_cell('stim_freq', paramnames);
    indstimd = find_ind_str_in_cell('stim_dur', paramnames);

    % Find indices of parameters that are influenced by other parameters
    indcpamp = find_ind_str_in_cell('cp_amp', paramnames);
    indcpper = find_ind_str_in_cell('cp_per', paramnames);
    indcpnum = find_ind_str_in_cell('cp_num', paramnames);

    % Update current pulse period (ms)
    paramvals(indcpper) = floor(1000/paramvals(indstimf));

    % Update number of current pulses                                        
    paramvals(indcpnum) = ceil(paramvals(indstimd)/paramvals(indcpper));
case 'noexp'
    % Do nothing
end

% Update current pulse amplitude (nA)
switch experimentname
case 'RTCl'
    paramvals(indcpamp) = 4*(paramvals(indREdia)/10)^2;
case 'm3ha'
    paramvals(indcpamp) = 0.2*(paramvals(indREdia)/10)^2;
case 'noexp'
    % Do nothing
end

%% Update GABAB-receptor-dependent parameters
switch experimentname
case 'm3ha'
    % Find indices of parameters that influence other parameters
    indpCond = find_ind_str_in_cell('pCond', paramnames);
    indgIncr = find_ind_str_in_cell('gIncr', paramnames);

    % Find indices of parameters that are influenced by other parameters
    indgabaaGmax = find_ind_str_in_cell('TCgabaaGmax', paramnames);
    indgababAmp = find_ind_str_in_cell('TCgababAmp', paramnames);
    indgababTrise = find_ind_str_in_cell('TCgababTrise', paramnames);
    indgababTfallFast = find_ind_str_in_cell('TCgababTfallFast', paramnames);
    indgababTfallSlow = find_ind_str_in_cell('TCgababTfallSlow', paramnames);
    indgababW = find_ind_str_in_cell('TCgababW', paramnames);
    
    % Update parameters that are influenced by other parameters
    pCond = paramvals(indpCond);
    gIncr = paramvals(indgIncr);
    if paramvals(indgabaaGmax) ~= 0             % Not in bicuculline mode
        paramvals(indgabaaGmax) = gtemp_gmax(pCond) * 2 * gIncr/100;  
                                    % update maximal GABA-A conductance (uS)
                                    %   2 times that of GABA-B 
                                    %   based on Huguenard & Prince, 1994
    end
    paramvals(indgababAmp) = gtemp_amp(pCond) * gIncr/100;
                                    % update GABA-B conductance amplitude (uS)
    paramvals(indgababTrise) = gtemp_Trise(pCond);
                                    % update rising phase time constant (ms)
    paramvals(indgababTfallFast) = gtemp_TfallFast(pCond);
                                    % update fast decay time constant (ms)
    paramvals(indgababTfallSlow) = gtemp_TfallSlow(pCond);
                                    % update slow decay time constant (ms)
    paramvals(indgababW) = gtemp_w(pCond);
                                    % update weight of the fast decay

case {'RTCl', 'noexp'}
    % Do nothing
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}
