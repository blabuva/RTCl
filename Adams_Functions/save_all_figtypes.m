function save_all_figtypes (fig, filename, varargin)
%% Save figures using all figure types provided
% Usage: save_all_figtypes (fig, filename, varargin)
% Arguments:
%       fig         - figure to save
%                   must be a a figure object or a Simulink block diagram
%       filename    - file name
%                   must be a character vector
%       figtypes    figure type(s) for saving; 
%                       e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%                   could be anything recognised by the built-in 
%                       saveas() function 
%                   (see isfigtype.m under Adams_Functions)
%                   default == 'png'
%
% Requires:
%       /home/Matlab/Adams_Functions/isfigtype.m
%
% Used by:    
%       /media/adamX/RTCl/raster_plot.m
%       /media/adamX/RTCl/single_neuron.m
%       /media/adamX/RTCl/tuning_curves.m
%       /media/adamX/RTCl/tuning_maps.m
%       /media/adamX/m3ha/optimizer4gabab/run_neuron_once_4compgabab.m
%       /home/Matlab/Adams_Functions/plot_tuning_curve.m
%       /home/Matlab/Adams_Functions/plot_tuning_map.m
%
% File History:
% 2017-05-09 Created by Adam Lu
% 2017-11-08 Replaced figbase with [figbase, '.', figtypes]
% 2018-05-08 Changed tabs to spaces and limited width to 80
% 

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

% Add required inputs to an Input Parser
addRequired(iP, 'fig');                         % figure to save
    % TODO: validation function %);
addRequired(iP, 'filename', ...                 % file name
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));

% Add optional inputs to the Input Parser
addOptional(iP, 'figtypes', 'png', ...          % figure type(s) for saving
    @(x) min(isfigtype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, fig, filename, varargin{:});
[~, figtypes] = isfigtype(iP.Results.figtypes, 'ValidateMode', true);

%% Save figure(s)
if ~isempty(figtypes)    % if at least one figtype is provided
    % Break down file name
    [directory, figbase, ~] = fileparts(filename);

    % Save as figtypes
    if iscell(figtypes)        % if many figtypes provided
        % Save figure as each figtype
        nfigtypes = numel(figtypes);    % number of figure types for saving
        for f = 1:nfigtypes
            saveas(fig, fullfile(directory, ...
                                [figbase, '.', figtypes{f}]), figtypes{f});
        end
    elseif ischar(figtypes)        % if only one figtype provided
        % Save figure as the figtype
        saveas(fig, fullfile(directory, [figbase, '.', figtypes]), figtypes);
    end
else            % if no figtype provided
    % Simply use saveas()
    saveas(fig, filename);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

