function plot_tuning_curve(pvalues, readout, cols_to_plot, pislog, plabel, readout_label, col_labels, xlimits, ylimits, figname, varargin)
%% Plot a 1-dimensional tuning curve
% USAGE: plot_tuning_curve(pvalues, readout, cols_to_plot, pislog, plabel, readout_label, col_labels, xlimits, ylimits, figname, varargin)
% Arguments: TODO: argument requirements
%       pvalues         - column vector of parameter values
%       readout         - matrix of readout values where each column is a readout vector
%       cols_to_plot    - (opt) columns of the readout matrix to plot
%                       default == 1
%       pislog          - (opt) whether parameter values are to be plotted log scaled
%                       default == 0
%       plabel          - (opt) label for the parameter changed, suppress by setting value to 'suppress'
%                       default == 'Parameter'
%       readout_label   - (opt) label for the readout, suppress by setting value to 'suppress'
%                       default == 'Readout'
%       col_labels      - (opt) labels for the columns in the readout matrix, suppress by setting value to {'suppress'}
%                       default == 'Column #1', 'Column #2', ...
%       xlimits         - (opt) axes limits for the parameter values, suppress by setting value to -1
%                       default == expand by a little bit
%       ylimits         - (opt) axes limits for the readout
%                       default == NIL
%       figname         - (opt) figure name for saving file
%                       default == NIL
%       varargin        - 'FigTypes': figure type(s) for saving; e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%                       could be anything recognised by the built-in saveas() function
%                       (see isfigtype.m under Adams_Functions)
%                       default == 'png'
%
% Requires:
%        /home/Matlab/Adams_Functions/isfigtype.m
%        /home/Matlab/Adams_Functions/save_all_figtypes.m
%
% Used by:
%        /media/adamX/RTCl/tuning_curves.m
%
% 2017-04-17 Moved from tuning_curves.m
% 2017-04-17 Simplified code
% 2017-04-17 Set default arguments
% 2017-04-17 Color map is now based on number of columns to plot
% 2017-05-09 Added 'FigTypes' as a parameter-value pair argument
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error('Not enough input arguments, type ''help plot_tuning_curve'' for usage');
end

% TODO: Add required inputs to an Input Parser
iP = inputParser;

% Use Input Parser for parameter-value pairs
addParameter(iP, 'SingleColor', [0, 0, 1]);        % Color when ncols_to_plot is 1
addParameter(iP, 'FigTypes', 'png', ...            % figure type(s) for saving; e.g., 'png', 'fig', or {'png', 'fig'}, etc.
    @(x) min(isfigtype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, varargin{:});
singlecolor = iP.Results.SingleColor;
[~, figtypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);

%% Extract number of columns
ncols = size(readout, 2);

%% Set defaults
if nargin < 3 || isempty(cols_to_plot)
    cols_to_plot = 1;
end
if nargin < 4 || isempty(pislog)
    pislog = 0;
end
if nargin < 5 || isempty(plabel)
    plabel = 'Parameter';
end
if nargin < 6 || isempty(readout_label)
    readout_label = 'Readout';
end
if nargin < 7 || isempty(col_labels)
    col_labels = cell(1, ncols);
    for c = 1:ncols
        col_labels{c} = ['Column #', num2str(c)];
    end
end

%% Create and clear figure if figname provided
if nargin >= 10 && ~isempty(figname)
    h = figure(floor(rand()*10^4)+1);
    clf(h);
end

%% Plot readout values against parameter values
numnoninf = sum(~isinf(readout), 1);    % number of parameter values that don't give infinite values
ncols_to_plot = length(cols_to_plot);   % number of columns to plot
cm = colormap(jet(ncols_to_plot));      % color map used
for c = 1:ncols_to_plot
    % Plot curve
    col = cols_to_plot(c);              % current column to plot
    if pislog
        % Note: can't have hold on before semilogx
        p = semilogx(pvalues, readout(:, col)); hold on;
    else
        p = plot(pvalues, readout(:, col)); hold on;
    end
    
    % Set color and display name
    if ncols_to_plot > 1
        set(p, 'Color', cm(c, :))
    elseif ncols_to_plot == 1
        set(p, 'Color', singlecolor);
    end
    if ~isequal(col_labels, {'suppress'})
        set(p, 'DisplayName', strrep(col_labels{col}, '_', '\_'));
    end

    % If there is only one value, mark with a circle
    if numnoninf(col) == 1
        set(p, 'Marker', 'o');
    end
end

%% Show legend only if readout has more than one columns
if ncols > 1
    legend('Location', 'eastoutside');
end

%% Restrict x axis if xlimits provided; otherwise expand the x axis by a little bit
if nargin >= 8 && ~isempty(xlimits)
    xlim(xlimits);
elseif ~isequal(xlimits, -1)
    xlim([pvalues(1) - (pvalues(2) - pvalues(1)), ...
        pvalues(end) + (pvalues(end) - pvalues(end-1))]);
end

%% Restrict y axis if ylimits provided
if nargin >= 9 && ~isempty(ylimits)
    ylim(ylimits);
end

%% Set title and axes labels
if ~isequal(plabel, 'suppress')
    xlabel(plabel);
end
if ~isequal(readout_label, 'suppress')
    ylabel(readout_label);
end
if ~isequal(plabel, 'suppress') && ~isequal(readout_label, 'suppress')
    title(strrep([readout_label, ' vs. ', plabel], '_', '\_'));
end

%% Save figure if figname provided
if nargin >= 10 && ~isempty(figname)
    save_all_figtypes(h, figname, figtypes);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

    if pislog
        p = semilogx(pvalues, readout, 'Color', 'b');
    else
        p = plot(pvalues, readout, 'Color', 'b');
    end
    if numnoninf == 1
        set(p, 'Marker', 'o');
    end

if ncols <= 6
    cm = colormap(prism);            % red, orange, yellow, green, blue, violet
elseif ncols <= 7
    cm = colormap(lines);            % blue, orange, yellow, purple, green, cyan, magenta
else
    cm = colormap(jet(ncols));        % color map used
end



%}
