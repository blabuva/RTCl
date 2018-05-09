function tuning_maps (infolder, outfolder, numactive, latency, actvel, nump, pnames, plabels, pislog, pvalues, nperp, ncells, actmode, loopmode, varargin)
%% TODO: Shows a tuning curve for numactive, latency, actvel, for each parameter changed
% USAGE: tuning_maps (infolder, outfolder, numactive, latency, actvel, nump, pnames, plabels, pislog, pvalues, nperp, ncells, actmode, loopmode, varargin)
% Arguments:
%    infolder    - the directory that contains the data to plot
%    outfolder    - (opt) the directory to place the output plots
%     TODO: other optional arguments
%    varargin    - 'FigTypes': figure type(s) for saving; e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%            could be anything recognised by the built-in saveas() function
%            (see isfigtype.m under Adams_Functions)
%            default == 'png'
%            - 'CellToPlot': cell to plot for the latency tuning map
%            default == ncells/2-1 when actmode ~= 3 and ncells/2-6 when actmode == 3
%
% Requires:
%        cd/raster_plot.m
%        infolder/*loopedparams.mat
%        /home/Matlab/Adams_Functions/isfigtype.m
%        /home/Matlab/Adams_Functions/save_all_figtypes.m
%        /home/Matlab/Adams_Functions/get_loopedparams.m
%        /home/Matlab/Adams_Functions/vec2array.m
%        /home/Matlab/Adams_Functions/plot_tuning_map.m
%
% Used by:
%        cd/raster_plot.m
%
% 2017-05-04 Adapted from tuning_curves.m
% 2017-05-09 Changed plotcurves to plottuning
% 2017-05-09 Changed 'FigType' to 'FigTypes', which now accepts a cell array of figure types
% 2017-05-09 Replaced saveas() with save_all_figtypes()
%

%% Set readout labels
numactive_label = 'Number of activated cells';
latency_label = 'Latency to activation (seconds)';
actvel_label = 'Activation velocity (cells/second)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add directories to search path for required functions
if exist(fullfile(pwd, '/Adams_Functions/'), 'dir') == 7
    functionsdirectory = pwd;
elseif exist('/home/Matlab/', 'dir') == 7
    functionsdirectory = '/home/Matlab/';
elseif exist('/scratch/al4ng/Matlab/', 'dir') == 7
    functionsdirectory = '/scratch/al4ng/Matlab/';
else
    error('Valid functionsdirectory does not exist!');
end
addpath(fullfile(functionsdirectory, '/Adams_Functions/'));        % for isfigtype.m, get_loopedparams.m,
                                    %     vec2array.m
                                    %    & plot_tuning_map.m
%% Check arguments
%%% TODO

%% Set outfolder and create if doesn't exist
if nargin < 2 || isempty(outfolder)
    outfolder = infolder;
end
if exist(outfolder, 'dir') ~= 7
    mkdir(outfolder);
end

%% Get numactive, latency, actvel if not provided
if nargin < 3 || isempty(numactive) || isempty(latency) || isempty(actvel)
    [~, numactive, latency, actvel] = ...
            raster_plot(infolder, 'OutFolder', outfolder, ...
                'RenewParpool', 0, 'PlotSpikes', 0, 'PlotTuning', 0);
end

%% Get looped parameters info if not provided
if nargin < 6 || isempty(nump) || isempty(pnames) || isempty(plabels) ...
        || isempty(pislog) || isempty(pvalues) || isempty(nperp) ...
        || isempty(ncells) || isempty(actmode) || isempty(loopmode)
    [nump, pnames, plabels, pislog, pvalues, nperp, ~, ~, ncells, actmode, loopmode] = get_loopedparams(infolder);
end

%% Use Input Parser for parameter-value pairs
iP = inputParser;
addParameter(iP, 'CellToPlot', -1);        % cell to plot for the latency tuning map
addParameter(iP, 'FigTypes', 'png', ...            % figure type(s) for saving; e.g., 'png', 'fig', or {'png', 'fig'}, etc.
    @(x) min(isfigtype(x, 'ValidateMode', true)));
parse(iP, varargin{:});
latency_cell_to_plot = iP.Results.CellToPlot;
[~, figtypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);

%% Find latency of just the cell to plot
if latency_cell_to_plot == -1
    switch actmode
        case {1, 2, 4, 5}
            latency_cell_to_plot = ncells/2 - 1;
        case 3
            latency_cell_to_plot = ncells/2 - 6;
    end
end
latency_this = latency(:, latency_cell_to_plot + 1);

%% Reorganize numactive, latency, actvel so that it is an array with nump dimensions
numactive = vec2array(numactive, nperp);
latency_this = vec2array(latency_this, nperp);
actvel = vec2array(actvel, nperp);

%% Create heat maps if there are two parameters
if nump == 2
    % Plot number of activated cells against parameter values
    figname = fullfile(outfolder, ['numactive_vs_', pnames{1}, '_', pnames{2}]);
    climits = [0, ncells];
    plot_tuning_map (pvalues, numactive, 'PisLog', pislog, 'PLabels', plabels, ...
            'ReadoutLabel', numactive_label, 'CLim', climits, ...
            'Figname', figname, 'FigTypes', figtypes)

    % Plot latency to activation (seconds) of cells against parameter values
    figname = fullfile(outfolder, ['latency_vs_', pnames{1}, '_', pnames{2}]);
    climits = [0, max(1, max(max(latency(~isinf(latency)))) + 1)];
    latency_label_this = [latency_label, ' for Cell#', num2str(latency_cell_to_plot)];
    plot_tuning_map (pvalues, latency_this, 'PisLog', pislog, 'PLabels', plabels, ...
            'ReadoutLabel', latency_label_this, 'CLim', climits, ...
            'Figname', figname, 'FigTypes', figtypes)


    % If not stimulating everywhere, plot activation velocity (cells/second) against parameter values
    if (actmode == 1 || actmode == 3 || actmode == 4)
        figname = fullfile(outfolder, ['actvel_vs_', pnames{1}, '_', pnames{2}]);
        climits = [0, max(max(actvel))+0.1];
        plot_tuning_map (pvalues, actvel, 'PisLog', pislog, 'PLabels', plabels, ...
                'ReadoutLabel', actvel_label, 'CLim', climits, ...
                'Figname', figname, 'FigTypes', figtypes)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}
