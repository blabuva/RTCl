function [pvalues, numactive, latency, actvel, nump, pnames, plabels, pislog, nperp, ncells, actmode, latency_cells_to_plot] = tuning_curves (infolders, outfolder, numactive, latency, actvel, nump, pnames, plabels, pislog, pvalues, nperp, ncells, actmode, loopmode, varargin)
%% Shows a tuning curve for numactive, latency, actvel, for each parameter changed
% USAGE: [pvalues, numactive, latency, actvel, nump, pnames, plabels, pislog, nperp, ncells, actmode, latency_cells_to_plot] = tuning_curves (infolders, outfolder, numactive, latency, actvel, nump, pnames, plabels, pislog, pvalues, nperp, ncells, actmode, loopmode, varargin)
% Arguments:
%    infolders    - the directory or a cell array of directories that contains the data to plot
%    outfolder    - (opt) the directory to place the output plots
%     TODO: other optional arguments; TODO: argument checks
%    varargin    - 'FigTypes': figure type(s) for saving; e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%            could be anything recognised by the built-in saveas() function
%            (see isfigtype.m under Adams_Functions)
%            default == 'png'
%            - 'CellsToPlot': cell(s) to plot for latencies
%            must be a nonnegative integer vector
%            default == 0:(ncells/2-1) if infolder is a cell array
%                    otherwise == ncells/2-1 when actmode ~= 3 and ncells/2-6 when actmode == 3
%            - 'PCondNames': names of parameters to be varied across conditions, e.g., {'stim_freq', 'tau'}
%            must be a cell array of strings or character arrays
%            default == {}
%    NOTE: If infolders is a cell array, 
%        'PCondNames' must be provided unless the same two parameters are varied inside each infolder
%
% Requires:
%        cd/raster_plot.m
%        infolder/*loopedparams.mat
%        /home/Matlab/Downloaded_Functions/dirr.m
%        /home/Matlab/Adams_Functions/isfigtype.m
%        /home/Matlab/Adams_Functions/save_all_figtypes.m
%        /home/Matlab/Adams_Functions/find_ind_str_in_cell.m
%        /home/Matlab/Adams_Functions/vec2cell.m
%        /home/Matlab/Adams_Functions/check_and_collapse_identical_contents.m
%        /home/Matlab/Adams_Functions/get_loopedparams.m
%        /home/Matlab/Adams_Functions/plot_tuning_curve.m
%
% Used by:
%        cd/raster_plot.m
%
% 2017-04-14 Moved from raster_plot.m
% 2017-04-17 Removed indices from the argument list of plot_tuning_curve()
% 2017-04-17 Moved plot_tuning_curve() to /home/Matlab/Adams_Functions/plot_tuning_curve.m
% 2017-04-17 Added the case when infolder is a cell array
% 2017-04-18 Added 'FigType' as a parameter-value pair argument
% 2017-04-18 Added 'CellsToPlot' as a parameter-value pair argument
% 2017-05-04 Added loopmode
% 2017-05-08 Used {'suppress'} to suppress column labels when not comparing across conditions
% 2017-05-09 Now pcondname can be a substring of paramname
% 2017-05-09 Changed plotcurves to plottuning
% 2017-05-09 Added figtypes and the function save_all_figtypes()
% 2017-05-09 Combined 'CellsToPlot' with 'CellToPlot'
% 2017-05-09 Combined 'FigType' with 'FigTypes'
% 2017-05-09 Replaced saveas() with save_all_figtypes()
% 2017-05-12 Changed reorganize_readout() to use vec2cell.m
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
addpath(fullfile(functionsdirectory, '/Downloaded_Functions/'));    % for dirr.m
addpath(fullfile(functionsdirectory, '/Adams_Functions/'));        % for isfigtype.m, save_all_figtypes.m,
                                    %    find_ind_str_in_cell.m, vec2cell.m
                                    %     check_and_collapse_identical_contents.m 
                                    %    & plot_tuning_curve.m & get_loopedparams.m

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error('Not enough input arguments, type ''help tuning_curves'' for usage');
end

% Add required inputs to an Input Parser
iP = inputParser;
%%% TODO

% Add optional inputs to the Input Parser
%%% TODO

%% Check to see if infolders is a cell array; if not, make it a cell array
if iscell(infolders)
    across_conditions = 1;        % whether the tuning curves across conditions are compared
    ncond = numel(infolders);    % number of conditions to compare
else
    infolders = {infolders};
    across_conditions = 0;
end

%% Set outfolder and create if doesn't exist
if nargin < 2 || isempty(outfolder)
    if across_conditions
        outfolder = 'tuning_curve_across_conditions';
    else
        outfolder = infolders{1};
    end
end
if exist(outfolder, 'dir') ~= 7
    mkdir(outfolder);
end

%% Get numactive, latency, actvel if not provided
if nargin < 3 || isempty(actvel)
    if across_conditions
        numactive = cell(1, ncond);
        latency = cell(1, ncond);
        actvel = cell(1, ncond);
        for k = 1:ncond
            [~, numactive{k}, latency{k}, actvel{k}] = ...
                    raster_plot(infolders{k}, 'OutFolder', outfolder, ...
                        'RenewParpool', 0, 'PlotSpikes', 0, 'PlotTuning', 0);
        end
    else
        [~, numactive, latency, actvel] = ...
                raster_plot(infolders{1}, 'OutFolder', outfolder, ...
                    'RenewParpool', 0, 'PlotSpikes', 0, 'PlotTuning', 0);
    end
end

%% Get looped parameters info if not provided
if nargin < 6 || isempty(actmode)        %%% TODO: might need more isempty()
    if across_conditions
        % Invariants for each condition
        nump = zeros(1, ncond);
        pnames = cell(1, ncond);
        plabels = cell(1, ncond);
        pislog = cell(1, ncond);
        ncells = zeros(1, ncond);
        actmode = zeros(1, ncond);    
        loopmode = cell(1, ncond);    

        % Variants for each condition
        pvalues = cell(1, ncond);    % parameter values might be different for each condition
        nperp = cell(1, ncond);        % number of parameter values might be different for each condition

        % Get the parameters from each condition
        for k = 1:ncond
            [nump(k), pnames{k}, plabels{k}, pislog{k}, pvalues{k}, nperp{k}, ...
                ~, ~, ncells(k), actmode(k), loopmode{k}] = get_loopedparams(infolders{k});
        end
    else
        [nump, pnames, plabels, pislog, pvalues, nperp, ~, ~, ncells, actmode, loopmode] = get_loopedparams(infolders{1});
    end
end

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'FigTypes', 'png', ...            % figure type(s) for saving; e.g., 'png', 'fig', or {'png', 'fig'}, etc.
    @(x) min(isfigtype(x, 'ValidateMode', true)));
addParameter(iP, 'CellsToPlot', 0:ncells/2-1, ...    % cells to plot for latencies when not comparing across conditions
    @(x) validateattributes(x, {'numeric'}, {'vector', 'nonnegative', 'integer'}));
addParameter(iP, 'PCondNames', {}, ...        % names of parameters to be varied across conditions
    @(x) assert(iscell(x) && (min(cellfun(@ischar, x)) || min(cellfun(@isstring, x))), ...
        'PCondNames must be a cell array of strings or character arrays!'));

% Read from the Input Parser
parse(iP, varargin{:});
[~, figtypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);
latency_cells_to_plot = iP.Results.CellsToPlot;
pcondnames = iP.Results.PCondNames;

%% Check if invariants are the same across conditions; if so, rename invariants
if across_conditions
    nump = check_and_collapse_identical_contents(nump, 'nump');
    pnames = check_and_collapse_identical_contents(pnames, 'pnames');
    plabels = check_and_collapse_identical_contents(plabels, 'plabels');
    pislog = check_and_collapse_identical_contents(pislog, 'pislog');
    ncells = check_and_collapse_identical_contents(ncells, 'ncells');
    actmode = check_and_collapse_identical_contents(actmode, 'actmode');
end

%% Create cell labels
cell_labels = cell(1, ncells);
for c = 1:ncells
    cell_labels{c} = ['Cell', num2str(c-1)];
end

%% Find varied parameters values for each condition
if across_conditions && isempty(pcondnames)
    pcondnames = cell(1, nump);
    if nump == 2        % currently only supports the case when 2 different parameters are looped and varied
        for p = 1:nump                % for each parameter that is varied within a condition
            % Find parameter that is varied across conditions (i.e., the other parameter that is looped)
            q = 3 - p;            % 2 -> 1; 1 -> 2 
            pcondnames{p} = pnames{q};    % the name of the parameter that is varied across conditions
        end
    else
        error('Please specify the names of the parameters varied across conditions!');
    end
end

%% Create condition labels
pcondvalues = cell(1, nump);
cond_labels = cell(1, nump);
cellcond_labels = cell(1, nump);
if across_conditions    
    for p = 1:nump                % for each parameter that is varied within a condition
        pcondvalues{p} = cell(1, ncond);
        cond_labels{p} = cell(1, ncond);
        cellcond_labels{p} = cell(1, ncond);
        for k = 1:ncond
            % Extract parameter values from sim_params file
            sim_params_files = dirr(fullfile(infolders{k}, ...
                ['sim_params_', pnames{p}, '*']), '.csv');    
                            % all sim_params files for the current looped parameter
            fid = fopen(fullfile(infolders{k}, sim_params_files(1).name));        % the first file suffices
            simfilecontent = textscan(fid, '%s %f %s', 'Delimiter', ',');
            paramnames = simfilecontent{1};
            params_val = simfilecontent{2};
            pcondvalues{p}{k} = params_val(find_ind_str_in_cell(pcondnames{p}, paramnames));

            % Create condition label
            cond_labels{p}{k} = [pcondnames{p}, ' = ', num2str(pcondvalues{p}{k})];

            % Create cell-condition labels
            cellcond_labels{p}{k} = cell(1, ncells);
            for c = 1:ncells
                cellcond_labels{p}{k}{c} = ...
                    ['Cell', num2str(c-1), ', ', ...
                    pcondnames{p}, ' = ', num2str(pcondvalues{p}{k})];
            end            
        end
    end
end

%% Reorganize numactive, latency, actvel so that each vector is for one parameter
if across_conditions
    for k = 1:ncond
        numactive{k} = reorganize_readout (numactive{k}, nump, nperp{k});
        latency{k} = reorganize_readout (latency{k}, nump, nperp{k});
        actvel{k} = reorganize_readout (actvel{k}, nump, nperp{k});
    end;
else
    numactive = reorganize_readout (numactive, nump, nperp);
    latency = reorganize_readout (latency, nump, nperp);
    actvel = reorganize_readout (actvel, nump, nperp);
end

%% Create tuning curves
for p = 1:nump
    % Define color map if many tuning curves are plotted
    if across_conditions
        cm = colormap(lines(ncond));    % color map used
    end

    % Plot number of activated cells against parameter values
    figname = fullfile(outfolder, ['numactive_vs_', pnames{p}]);
    ylimits = [0, ncells];
    if across_conditions
        h = figure(floor(rand()*10^4)+1);
        clf(h);
        for k = 1:ncond
            plot_tuning_curve(pvalues{k}{p}, numactive{k}{p}, 1, ...
                        pislog(p), plabels{p}, numactive_label, {cond_labels{p}{k}}, ...
                        [], ylimits, [], 'SingleColor', cm(k, :));        
        end
        legend('Location', 'northeast');
        save_all_figtypes(h, figname, figtypes);
    else
        plot_tuning_curve(pvalues{p}, numactive{p}, 1, ...
                    pislog(p), plabels{p}, numactive_label, ...
                    {'suppress'}, [], ylimits, figname, ...
                    'FigTypes', figtypes);
        %% TODO: Pass 'figtypes' to plot_tuning_curve.m
    end

    % Plot latency to activation (seconds) of cells against parameter values
    figname = fullfile(outfolder, ['latency_vs_', pnames{p}]);
    if across_conditions
        h = figure(floor(rand()*10^4)+1);
        clf(h);
        for k = 1:ncond
            if latency_cells_to_plot == -1
                switch actmode
                    case {1, 2, 4, 5}
                        cols_to_plot = ncells/2 - 1 + 1;
                    case 3
                        cols_to_plot = ncells/2 - 6 + 1;
                end
            else
                cols_to_plot = latency_cells_to_plot + 1;
            end
            plot_tuning_curve(pvalues{k}{p}, latency{k}{p}, cols_to_plot, ...
                        pislog(p), plabels{p}, latency_label, cellcond_labels{p}{k}, ...
                        [], [], [], 'SingleColor', cm(k, :));
        end
        legend('Location', 'northeast');
        save_all_figtypes(h, figname, figtypes)
    else
        if latency_cells_to_plot == -1
            cols_to_plot = (0:(ncells/2-1)) + 1;
        else
            cols_to_plot = latency_cells_to_plot + 1;
        end
        ylimits = [0, max(1, max(latency{p}(~isinf(latency{p}))) + 1)];
        plot_tuning_curve(pvalues{p}, latency{p}, cols_to_plot, ...
                    pislog(p), plabels{p}, latency_label, ...
                    cell_labels, [], ylimits, figname, ...
                    'FigTypes', figtypes);
    end

    % If not stimulating everywhere, plot activation velocity (cells/second) against parameter values
    if (actmode == 1 || actmode == 3 || actmode == 4)
        figname = fullfile(outfolder, ['actvel_vs_', pnames{p}]);
        if across_conditions
            h = figure(floor(rand()*10^4)+1);
            clf(h);
            for k = 1:ncond
                plot_tuning_curve(pvalues{k}{p}, actvel{k}{p}, 1, ...
                        pislog(p), plabels{p}, actvel_label, {cond_labels{p}{k}}, ...
                        [], [], [], 'SingleColor', cm(k, :));
            end
            legend('Location', 'northeast');
            save_all_figtypes(h, figname, figtypes)
        else
            ylimits = [0, max(actvel{p})+0.1];
            plot_tuning_curve(pvalues{p}, actvel{p}, 1, ...
                    pislog(p), plabels{p}, actvel_label, ...
                    {'suppress'}, [], ylimits, figname, ...
                    'FigTypes', figtypes);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function readout_reorg = reorganize_readout (readout, nump, nperp)
%% Reorganize readout vector or matrix

%% Create class vector for the parameter classes
classvector = [];
for p = 1:nump
    classvector = [classvector; p * ones(nperp(p), 1)];
end

%% Reorganize readout vector into a cell array of vectors according to parameter class
readout_reorg = vec2cell(readout, classvector);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}
