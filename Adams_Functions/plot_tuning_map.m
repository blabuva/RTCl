function plot_tuning_map (pvalues, readout, varargin)
%% Plot a 2-dimensional tuning map
% Usage: plot_tuning_map (pvalues, readout, varargin)
% Arguments:
%       pvalues     - 2 vectors of parameter values
%                   must be a 2-element cell array of numeric vectors
%       readout     - a readout matrix corresponding to the parameter-pairs
%                   must be a 2-dimensional numeric array with dimensions 
%                   given by the lengths of elements of pvalues
%       varargin    - 'PisLog': whether parameter values are to be plotted log scaled
%                   must be a 2-element binary array
%                   default == [false, false];
%                   - 'PLabels': labels for parameters
%                   must be a 2-element cell array of string scalars or character vectors
%                   default == {'Parameter1', 'Parameter2'}
%                   - 'ReadoutLabel': label for readout matrix
%                   must be a string scalar or a character vector
%                   default == 'Readout'
%                   - 'XLim': limits of x axis
%                   must be a 2-element increasing numeric vector
%                   default == []
%                   - 'YLim': limits of y axis
%                   must be a 2-element increasing numeric vector
%                   default == []
%                   - 'CLim': limits of color axis
%                   must be a 2-element increasing numeric vector
%                   default == []
%                   - 'FigName': figure name for saving
%                   must be a string scalar or a character vector
%                   default == ''
%                   - 'FigTypes': figure type(s) for saving; e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%                   could be anything recognised by the built-in saveas() function
%                   (see isfigtype.m under Adams_Functions)
%                   default == 'png'
%
% Requires:
%       /home/Matlab/Adams_Functions/isfigtype.m
%       /home/Matlab/Adams_Functions/save_all_figtypes.m
%
% Used by:    
%       /media/adamX/RTCl/tuning_maps.m
%
% File History:
% 2017-05-05 Created by Adam Lu
% 2017-05-09 Added 'FigTypes' as a parameter-value pair argument
% 2017-12-18 Changed tabs to spaces
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error('Not enough input arguments, type ''help plot_tuning_map'' for usage');
end

% Add required inputs to an Input Parser
iP = inputParser;
addRequired(iP, 'pvalues', ...              % 2 vectors of parameter values
    @(x) assert(iscell(x) && numel(x) == 2 ...
        && isnumeric(x{1}) && isvector(x{1}) ...
        && isnumeric(x{2}) && isvector(x{1}), ...
        'pvalues must be a 2-element cell array of numeric vectors!'));
addRequired(iP, 'readout', ...              % a readout matrix corresponding to the parameter-pairs
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'PisLog', [0, 0], ...      % whether parameter values are to be plotted log scaled
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary', 'vector', 'numel', 2}));
addParameter(iP, 'PLabels', {'Parameter1', 'Parameter2'}, ...    % labels for parameters
    @(x) assert(iscell(x) && numel(x) == 2 ...
        && (min(cellfun(@ischar, x)) || min(cellfun(@isstring, x))), ...
        'pvalues must be a 2-element cell array of string scalars or character vectors!'));
addParameter(iP, 'ReadoutLabel', 'Readout', ...     % label for readout matrix
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
addParameter(iP, 'XLim', [], ...            % limits of x axis
    @(x) validateattributes(x, {'numeric'}, {'increasing', 'vector', 'numel', 2}));
addParameter(iP, 'YLim', [], ...            % limits of y axis
    @(x) validateattributes(x, {'numeric'}, {'increasing', 'vector', 'numel', 2}));
addParameter(iP, 'CLim', [], ...            % limits of color axis
    @(x) validateattributes(x, {'numeric'}, {'increasing', 'vector', 'numel', 2}));
addParameter(iP, 'FigName', '', ...         % figure name for saving
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
addParameter(iP, 'FigTypes', 'png', ...     % figure type(s) for saving
    @(x) min(isfigtype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, pvalues, readout, varargin{:});
pislog = iP.Results.PisLog;
plabels = iP.Results.PLabels;
readout_label = iP.Results.ReadoutLabel;
xlimits = iP.Results.XLim;
ylimits = iP.Results.YLim;
climits = iP.Results.CLim;
figname = iP.Results.FigName;
[~, figtypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);

% Check relationships between arguments
if size(readout) ~= cellfun(@length, pvalues)
    error('Readout matrix must have dimensions given by the lengths of elements of pvalues!');
end

%% Create and clear figure if figname provided
if ~isempty(figname)
    h = figure(floor(rand()*10^4)+1);
    clf(h);
end

%% Plot readout values against parameter values
surf(pvalues{1}, pvalues{2}, readout');
if pislog(1)
    set(gca, 'XScale', 'log');
elseif pislog(2)
    set(gca, 'YScale', 'log');
end
set(gca, 'CLim', climits);
colorbar;

%% Restrict x axis if xlimits provided
if ~isempty(xlimits)
    xlim(xlimits);
end

%% Restrict y axis if ylimits provided
if ~isempty(ylimits)
    ylim(ylimits);
end

%% Set title and axes labels; allow to be suppressed
if ~isequal(plabels{1}, 'suppress')
    xlabel(plabels{1});
end
if ~isequal(plabels{2}, 'suppress')
    ylabel(plabels{2});
end
if ~isequal(plabels{1}, 'suppress') && ~isequal(plabels{2}, 'suppress') && ~isequal(readout_label, 'suppress')
    title(strrep(readout_label, '_', '\_'));
%    title(strrep([readout_label, ' vs. ', plabels{1}, ' and ', plabels{2}], '_', '\_'));
end

%% Save figure if figname provided
if ~isempty(figname)
    save_all_figtypes(h, figname, figtypes);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% imagesc([pvalues{1}(1), pvalues{1}(end)], [pvalues{2}(1), pvalues{2}(end)], readout);
% imagesc([pvalues{1}(1), pvalues{1}(end)], [pvalues{2}(end), pvalues{2}(1)], flipud(readout));

    saveas(h, figname);

%}
