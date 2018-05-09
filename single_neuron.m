function [data] = single_neuron (infolder, varargin)
%% Shows single neuron traces for different neurons or for different properties in the same neuron
% USAGE: [data] = single_neuron (infolder, varargin)
% Arguments:
%       infolder    - the name of the directory containing the .syn files, e.g. '20170317T1127_Ggaba_0.01'
%                    must be a character array
%       varargin    - 'FigTypes': figure type(s) for saving; e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%                   could be anything recognised by the built-in saveas() function
%                   (see isfigtype.m under Adams_Functions)
%                   default == 'png'
%                   - 'OutFolder': the name of the directory that the plots will be placed
%                   must be a directory
%                   default: same as infolder
%                   - 'MaxNumWorkers': maximum number of workers for running NEURON 
%                   must a positive integer
%                   default: 20
%                   - 'RenewParpool': whether to renew parpool every batch to release memory
%                   must be logical 1 (true) or 0 (false)
%                   default: true, but alway false if plotspikes == false
%                   - 'CellsToPlot': the ID #s for cells whose voltage & chloride concentration traces are to be plotted
%                   must be a numeric array with elements that are integers between 0 and ncells
%                   default: [act, act_left1, act_left2, far], whose values are saved in sim_params.csv
%                   - 'PropertiesToPlot': property #s of special neuron to record to be plotted 
%                   maximum range: 1~12, must be consistent with net.hoc
%                   must be a numeric array with elements that are integers between 0 and 12
%                   default: 1:12
%                   legend: 1 - voltage (mV) trace
%                           2 - sodium current (mA/cm2) trace
%                           3 - potassium current (mA/cm2) trace
%                           4 - calcium current (mA/cm2) trace 
%                           5 - calcium concentration (mM) trace
%                           6 - GABA-A chloride current (nA) trace
%                           7 - GABA-A bicarbonate current (nA) trace
%                           8 - chloride current (mA/cm2) trace
%                           9 - chloride concentration (mM) trace
%                           10 - chloride concentration (mM) in inner annuli trace
%                           11 - chloride reversal potential trace
%                           12 - GABA-A reversal potential trace
%                   NOTE: must be consistent with proplabels & net.hoc
%
% Requires:
%        infolder/*.singv OR infolder/*.singcli OR infolder/*.singsp
%        infolder/['sim_params_', pstring, '.csv'] 
%           for all the possible parameter strings
%        /home/Matlab/Downloaded_Functions/dirr.m
%        /home/Matlab/Downloaded_Functions/subaxis.m
%        /home/Matlab/Adams_Functions/isfigtype.m
%        /home/Matlab/Adams_Functions/save_all_figtypes.m
%        /home/Matlab/Adams_Functions/find_ind_str_in_cell.m
% Used by:
%        cd/neuronlaunch.m
%
% 2017-02-13 Modified from raster_plot.m
% 2017-02-23 Now plots voltage traces of four neurons at once
% 2017-02-26 Now plots chloride concentration traces too
% 2017-03-04 Added traces for the special neuron
% 2017-03-05 Added filetype
% 2017-03-07 Changed v_figname to v_figsuffix, etc.
% 2017-03-07 Now deals with 2 special neurons for each trial
% 2017-03-13 Added line for stimulation
% 2017-03-13 Mark stimulation with dotted lines instead
% 2017-03-14 labels should now include the full length
% 2017-03-16 Time unit labels are now flexible
% 2017-03-17 Changed arguments: removed filetype and added IDs, IDlabels, properties, proplabels as arguments
% 2017-03-17 Now makes parfor more effective by separating all plots into different iterations
% 2017-03-17 Now plots voltage & chloride concentration heat maps for each set of neurons
% 2017-03-27 Filenames are now modified in neuronlaunch.m
% 2017-03-27 figsuffix now contains .png
% 2017-03-30 Heat maps now have the same colorbar limits
% 2017-03-30 Spike data are now overlayed on heat maps
% 2017-03-31 Now reads simulation parameters from the csv file
% 2017-04-04 Moved IDlabels and proplabels from neuronlaunch.m
% 2017-04-04 REuseca is now an argument of plot_single_neuron_data()
% 2017-04-04 Removed IDlabels and proplabels and create them inside plot_single_neuron_data()
% 2017-04-06 Added close all inside parfor
% 2017-04-07 Added zoom3.png
% 2017-04-10 Added renewparpool
% 2017-04-10 Made renewparpool, maxnumworkers arguments
% 2017-04-18 Added 'FigType' as a parameter-value pair argument
% 2017-05-09 Changed 'FigType' to 'FigTypes', which now accepts a cell array of figure types
% 2017-05-09 Replaced saveas() with save_all_figtypes()
% 2017-05-09 Revamped Input Parser scheme
%

%% Set parameters
nzooms = 4;            % number of different time intervals to plot for each data
nppf = nzooms + 1;        % number of plots per file

%% Set figure name suffices
v_figsuffix = '_selected_soma_voltage.png';
cli_figsuffix = '_selected_soma_cli.png';
sp_figsuffix = '_alltraces.png';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error('An infolder is required, type ''help single_neuron'' for usage');
end

% Add required inputs to an Input Parser
iP = inputParser;
addRequired(iP, 'infolder', @isdir);            % the name of the directory containing the .syn files

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'FigTypes', 'png', ...            % figure type(s) for saving; e.g., 'png', 'fig', or {'png', 'fig'}, etc.
    @(x) min(isfigtype(x, 'ValidateMode', true)));
addParameter(iP, 'OutFolder', '@infolder', @isdir);    % the name of the directory that the plots will be placed
addParameter(iP, 'MaxNumWorkers', 20, ...        % maximum number of workers for running NEURON
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'RenewParpool', true, ...        % whether to renew parpool every batch to release memory
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'CellsToPlot', [], ...            % the ID #s for cells to be plotted
    @(x) validateattributes(x, {'numeric'}, {'vector', 'nonnegative', 'integer', '>=', 0}));
addParameter(iP, 'PropertiesToPlot', 1:1:12, ...    % property #s of special neuron to record to be plotted
    @(x) validateattributes(x, {'numeric'}, {'vector', 'positive', 'integer', ...
                            '>', 0, '<', 13}));

% Read from the Input Parser
parse(iP, infolder, varargin{:});
[~, figtypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);
outfolder = iP.Results.OutFolder;
maxnumworkers = iP.Results.MaxNumWorkers;
renewparpool = iP.Results.RenewParpool;
cellstoplot = iP.Results.CellsToPlot;
propertiestoplot = iP.Results.PropertiesToPlot;

% Change default arguments if necessary
if strcmp(outfolder, '@infolder')
    outfolder = infolder;
end

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
addpath(fullfile(functionsdirectory, '/Downloaded_Functions/'));    % for dirr.m & subaxis.m
addpath(fullfile(functionsdirectory, '/Adams_Functions/'));        % for isfigtype.m & find_ind_str_in_cell.m

%% Find all .singv, .singcli or .singsp files
files = dirr(infolder, '.sing*');
nfiles = length(files);

%% Set up plots
nplots = nfiles * nppf;
filenames = cell(1, nplots);        % stores file name for each plot
filetypes = cell(1, nplots);        % stores file type for each plot
fullfignames = cell(1, nplots);        % stores full figure names for each plot
ncells = zeros(1, nplots);        % stores ncells for each plot
tstarts = zeros(1, nplots);        % stores tstart for each plot
tstops = zeros(1, nplots);        % stores tstop for each plot
stim_starts = zeros(1, nplots);        % stores stim_start for each plot
stim_durs = zeros(1, nplots);        % stores stim_dur for each plot
stim_freqs = zeros(1, nplots);        % stores stim_freq for each plot
REuseca = zeros(1, nplots);        % stores REuseca for each plot
act = zeros(1, nplots);            % stores act for each plot
act_left1 = zeros(1, nplots);        % stores act_left1 for each plot
act_left2 = zeros(1, nplots);        % stores act_left2 for each plot
far = zeros(1, nplots);            % stores far for each plot
IDs_cell = cell(1, nplots);        % stores default neuron ID #s to plot for each plot
for i = 1:nfiles
    % Set things common for all time intervals to plot
    for j = 1:nppf
        % Find current index of plots
        ci = nppf*(i-1)+j;            % current index of plots

        % Store file name
        filenames{ci} = files(i).name;

        % Set filetype according to filename
        if strfind(files(i).name, '.singv')
            filetypes{ci} = 'v';
        elseif strfind(files(i).name, '.singcli')
            filetypes{ci} = 'cli';
        elseif strfind(files(i).name, '.singsp')
            filetypes{ci} = 'sp';
        end

        % Extract parameters from sim_params file
        [~, filebase, ~] = fileparts(files(i).name);
        temparray = strsplit(filebase, '_');
        simfilename = ['sim_params_', strjoin(temparray(2:end), '_'), '.csv'];
        fid = fopen(fullfile(infolder, simfilename));
        simfilecontent = textscan(fid, '%s %f %s', 'Delimiter', ',');
        paramnames = simfilecontent{1};
        params_val = simfilecontent{2};
        ncells(ci) = params_val(find_ind_str_in_cell('ncells', paramnames, 'SearchMode', 'exact'));
        tstarts(ci) = params_val(find_ind_str_in_cell('tstart', paramnames, 'SearchMode', 'exact'));
        tstops(ci) = params_val(find_ind_str_in_cell('tstop', paramnames, 'SearchMode', 'exact'));
        stim_starts(ci) = params_val(find_ind_str_in_cell('stim_start', paramnames, 'SearchMode', 'exact'));
        stim_durs(ci) = params_val(find_ind_str_in_cell('stim_dur', paramnames, 'SearchMode', 'exact'));
        stim_freqs(ci) = params_val(find_ind_str_in_cell('stim_freq', paramnames, 'SearchMode', 'exact'));
        REuseca(ci) = params_val(find_ind_str_in_cell('REuseca', paramnames, 'SearchMode', 'exact'));
        act(ci) = params_val(find_ind_str_in_cell('act', paramnames, 'SearchMode', 'exact'));
        act_left1(ci) = params_val(find_ind_str_in_cell('act_left1', paramnames, 'SearchMode', 'exact'));
        act_left2(ci) = params_val(find_ind_str_in_cell('act_left2', paramnames, 'SearchMode', 'exact'));
        far(ci) = params_val(find_ind_str_in_cell('far', paramnames, 'SearchMode', 'exact'));
        fclose(fid);

        if ~isempty(cellstoplot)
            IDs_cell{ci} = cellstoplot;
        else
            IDs_cell{ci} = [act(ci), act_left1(ci), ...
                act_left2(ci), far(ci)];    % default ID #s for neurons whose voltage is to be plotted
        end

    end

    % Set general figure names according to filename
    if strfind(files(i).name, '.singv')
        figname = strrep(files(i).name, '.singv', v_figsuffix);
    elseif strfind(files(i).name, '.singcli')
        figname = strrep(files(i).name, '.singcli', cli_figsuffix);    
    elseif strfind(files(i).name, '.singsp')
        figname = strrep(files(i).name, '.singsp', sp_figsuffix);
    end

    % Create full figure names with modifications
    fullfignames{nppf*(i-1)+1} = fullfile(outfolder, figname);
    fullfignames{nppf*(i-1)+2} = strrep(fullfignames{nppf*(i-1)+1}, '.png', '_zoom1.png');
    fullfignames{nppf*(i-1)+3} = strrep(fullfignames{nppf*(i-1)+1}, '.png', '_zoom2.png');
    fullfignames{nppf*(i-1)+4} = strrep(fullfignames{nppf*(i-1)+1}, '.png', '_zoom3.png');
    fullfignames{nppf*(i-1)+5} = strrep(fullfignames{nppf*(i-1)+1}, 'selected', 'heatmap');

    % Create time limits for different time intervals to plot
    tstarts(nppf*(i-1)+1) = tstarts(ci);                        % 0 ms
    tstops(nppf*(i-1)+1) = tstops(ci);                        % 233000 ms
    tstarts(nppf*(i-1)+2) = max(stim_starts(ci)*2/3, tstarts(ci));            % 2000 ms
%    tstops(nppf*(i-1)+2) = min(max(tstops(ci)*10/63, stim_starts(ci)*5/3), tstops(ci));    % 10000 ms
    tstops(nppf*(i-1)+2) = min(max(10000, stim_starts(ci)*5/3), tstops(ci));    % 10000 ms
    tstarts(nppf*(i-1)+3) = max(stim_starts(ci)*29/30, tstarts(ci));            % 2900 ms
%    tstops(nppf*(i-1)+3) = min(max(tstops(ci)*4/63, stim_starts(ci)*32/30), tstops(ci));    % 4000 ms
    tstops(nppf*(i-1)+3) = min(max(4000, stim_starts(ci)*32/30), tstops(ci));    % 4000 ms
    tstarts(nppf*(i-1)+4) = max(stim_starts(ci) + stim_durs(ci) - 1000, tstarts(ci));    % 202000 ms
    tstops(nppf*(i-1)+4) = min(stim_starts(ci) + stim_durs(ci) + 100, tstops(ci));        % 203300 ms
    tstarts(nppf*(i-1)+5) = tstarts(ci);                        % 0 ms
    tstops(nppf*(i-1)+5) = tstops(ci);                        % 233000 ms
end

%% Create plots
data = cell(1, nplots);                % some elements will be empty but the indexing is necessary for parfor
ct = 0;                        % counts number of trials completed
poolobj = gcp('nocreate');            % get current parallel pool object without creating a new one
if isempty(poolobj)
    poolobj = parpool;            % create a default parallel pool object
    oldnumworkers = poolobj.NumWorkers;    % number of workers in the default parallel pool object
else
    oldnumworkers = poolobj.NumWorkers;    % number of workers in the current parallel pool object
end
numworkers = min(oldnumworkers, maxnumworkers);    % number of workers to use for running NEURON
if renewparpool
    delete(poolobj);            % delete the parallel pool object to release memory
end
while ct < nplots            % while not trials are completed yet
    first = ct + 1;            % first trial in this batch
    if renewparpool && ct + numworkers <= nplots    % if memory is to be released
        last = ct + numworkers;            % limit the batch to numworkers
    else
        last = nplots;
    end
    if renewparpool
        poolobj = parpool('local', numworkers);    % recreate a parallel pool object 
                            % using fewer workers to prevent running out of memory
    end
    parfor k = first:last
        if strcmp(filetypes{k}, 'v') || strcmp(filetypes{k}, 'cli')
            % Plot voltage or chloride concentration traces
            if mod(k, nppf) == 0
                plot_heat_map(filetypes{k}, tstarts(k), tstops(k), filenames{k}, fullfignames{k}, ...
                        stim_starts(k), stim_durs(k), stim_freqs(k), infolder, figtypes);
            elseif mod(k, nppf) == 1
    %            data{k} = ...
                    plot_single_neuron_data(filetypes{k}, ncells(k), REuseca(k), ...
                        tstarts(k), tstops(k), filenames{k}, fullfignames{k}, ...
                        IDs_cell{k}, stim_starts(k), stim_durs(k), stim_freqs(k), infolder, figtypes);
            else
                plot_single_neuron_data(filetypes{k}, ncells(k), REuseca(k), ...
                    tstarts(k), tstops(k), filenames{k}, fullfignames{k}, ...
                    IDs_cell{k}, stim_starts(k), stim_durs(k), stim_freqs(k), infolder, figtypes);
            end
        elseif strcmp(filetypes{k}, 'sp')
            % Plot other properties traces for special neurons
            if mod(k, nppf) == 0
                % No heat map; do nothing
            elseif mod(k, nppf) == 1
    %            data{k} = ...
                    plot_single_neuron_data(filetypes{k}, ncells(k), REuseca(k), ...
                        tstarts(k), tstops(k), filenames{k}, fullfignames{k}, ...
                        propertiestoplot, stim_starts(k), stim_durs(k), stim_freqs(k), ...
                        infolder, figtypes);
            else
                plot_single_neuron_data(filetypes{k}, ncells(k), REuseca(k), ...
                    tstarts(k), tstops(k), filenames{k}, fullfignames{k}, ...
                    propertiestoplot, stim_starts(k), stim_durs(k), stim_freqs(k), infolder, figtypes);
            end
        end

        close all;
    end
    if renewparpool
        delete(poolobj);            % delete the parallel pool object to release memory
    end
    ct = last;                % update number of trials completed
end
if renewparpool
    poolobj = parpool('local', oldnumworkers);    % recreate a parallel pool object using the previous number of workers
end

%% Remove empty elements from data
data = data(~cellfun(@isempty, data));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = plot_heat_map(filetype, tstart, tstop, filename, figname, stim_start, stim_dur, stim_freq, infolder, figtypes)
%% Plot heat map

% Load single neuron data
data = load(fullfile(infolder, filename));

% Check (ID #s of neurons) to plot
ncells = size(data, 2) - 1;        % total number of columns in the data minus the time vector
if ncells < 1
    fprintf('Warning: No neurons to plot for this file!\n');
    return;
end

% Find range of data values to plot and spike time data filename
if strcmp(filetype, 'v')
    climits = [-100, 50];
    spifilename = strrep(filename, '.singv', '.spi');
elseif strcmp(filetype, 'cli')
    climits = [0, 100];
    spifilename = strrep(filename, '.singcli', '.spi');
end

% Load spike data
spikes = load(fullfile(infolder, spifilename));
if ~isempty(spikes)
    spikecelln = spikes(:, 1);
end

% Change units of time axis from ms to s if the total time is > 10 seconds
if ~isempty(spikes)
    [timevec, timelabel, xlim1, xlim2, stim_start_plot, stim_dur_plot, spiketimes] ...
        = set_time_units(data, tstart, tstop, stim_start, stim_dur, spikes);
else
    [timevec, timelabel, xlim1, xlim2, stim_start_plot, stim_dur_plot, ~] ...
        = set_time_units(data, tstart, tstop, stim_start, stim_dur);
end

% Find maximum and minimum time to plot
xmin = min([xlim1, min(timevec)]);
xmax = max([xlim2, max(timevec)]);

% Create plot
%h = figure(floor(rand()*10^4+, 1));
h = figure(10000);
clf(h);
hold on;

% Create heat map
imagesc([min(timevec), max(timevec)], [ncells-1, 0], flipud(data(:, 2:end)'));
set(gca, 'CLim', climits);
%HeatMap(flipud(data(:, 2:end)'), 'ColumnLabels', timevec, 'RowLabels', 0:ncells-1);        % Doesn't seem to work
%heatmap(timevec, 0:ncells-1, flipud(data(:, 2:end)'));        % Not available until R2017a
if ~isempty(spikes)
    plot(spiketimes, spikecelln, 'r.', 'MarkerSize', 1);                % plot spikes
end
line([stim_start_plot, stim_start_plot], [-1, ncells], ...
    'Color', 'r', 'LineStyle', '--');            % line for stimulation on
text(stim_start_plot + 0.5, ncells*0.95, ['Stim ON: ', num2str(stim_freq), ' Hz'], ...
    'Color', 'r');                        % text for stimulation on
line([stim_start_plot + stim_dur_plot, stim_start_plot + stim_dur_plot], ...
    [-1,  ncells], 'Color', 'r', 'LineStyle', '--');    % line for stimulation off
text(stim_start_plot + stim_dur_plot + 0.5, ncells*0.95, 'Stim OFF', ...
    'Color', 'r');                        % text for stimulation off
xlim([xmin, xmax]);                        % time range to plot
ylim([-1, ncells]);                        % cell ID runs from 0 to ncells-1
xlabel(timelabel);
ylabel('Neuron number');
colorbar;
if strcmp(filetype, 'v')
    title(['Somatic voltage (mV) for ', strrep(filename, '_', '\_')]);
elseif strcmp(filetype, 'cli')
    title(['Chloride concentration (mM) for ', strrep(filename, '_', '\_')]);
end

% Save figure
save_all_figtypes(h, figname, figtypes);
%close(h);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = plot_single_neuron_data(filetype, ncells, REuseca, tstart, tstop, filename, figname, ToPl, stim_start, stim_dur, stim_freq, infolder, figtypes)
%% Plot single neuron data

%% Extract info from arguments
nsubplots = length(ToPl);                % number of column vectors to plot

%% Create legend labels
if strcmp(filetype, 'v') || strcmp(filetype, 'cli')
    % Check ID numbers
    if max(ToPl) > ncells
        error('IDs are out of range!');
    end

    % Create ID labels
    labels = cell(1, ncells);    % ID labels
    for id = 0:ncells-1
        labels{id+1} = sprintf('Cell #%d', id);
    end
elseif strcmp(filetype, 'sp')
    % Set labels for properties of the special neuron to be plotted, must be consistent with net.hoc
    labels = { ...
        'v (mV)', 'ina (mA/cm2)', 'ik (mA/cm2)', ...
        'ica (mA/cm2)', 'cai (mM)', ...
        'Gicl (nA)', 'Gihco3 (nA)', ...
        'icl (mA/cm2)', 'cli (mM)', 'cli1 (mM)', ...
        'ecl (mV)', 'eGABA (mV)'};
    if REuseca == 0
        labels{find_ind_str_in_cell('ica (mA/cm2)', labels)} = 'drive_ch (um2 mM/ms)';
        labels{find_ind_str_in_cell('cai (mM)', labels)} = 'drive_ex (um2 mM/ms)';
    end
end


% Load data
data = load(fullfile(infolder, filename));

% Check (ID #s of neurons) or (properties of special neuron) to plot
ncols = size(data, 2) - 1;        % total number of columns in the data minus the time vector
if ncols < 1
    fprintf('Warning: No neurons or properties to plot for this file!\n');
    return;
end
for k = 1:nsubplots
    if strcmp(filetype, 'v') || strcmp(filetype, 'cli')
        if ToPl(k) < 0 || ToPl(k) > ncols - 1
            fprintf('Warning: ToPl(%d) is out of range; plotting first neuron or property instead\n', k);
            ToPl(k) = 0;
        end
    elseif strcmp(filetype, 'sp')
        if ToPl(k) < 1 || ToPl(k) > ncols
            fprintf('Warning: ToPl(%d) is out of range; plotting first neuron or property instead\n', k);
            ToPl(k) = 0;
        end
    end
end

% Create figure
%h = figure(floor(rand()*10^4+, 1));
h = figure(10000);
clf(h);

% Plot voltage trace for each neuron with iD # in ToPl
for k = 1:nsubplots
    % Change units of time axis from ms to s if the total time is > 10 seconds
    [timevec, timelabel, xlim1, xlim2, stim_start_plot, stim_dur_plot] = set_time_units(data, tstart, tstop, stim_start, stim_dur);

    % Generate subplot
    subaxis(nsubplots, 1, k, 'SpacingVert', 0.015)
    hold on;
    if strcmp(filetype, 'v') || strcmp(filetype, 'cli')
        plot(timevec, data(:, ToPl(k)+2), 'b', ...            % trace for neuron #i is in the i+2nd column
            'DisplayName', strrep(labels{ToPl(k)+1}, '_', '\_'));    % label for neuron #i is in the i+1st entry
    elseif strcmp(filetype, 'sp')
        plot(timevec, data(:, ToPl(k)+1), 'b', ...            % property #i is in the i+1st column
            'DisplayName', strrep(labels{ToPl(k)}, '_', '\_'));
    end
    xlim([xlim1, xlim2]);
    if strcmp(filetype, 'v')
        ylim([-100, 60]);
    elseif strcmp(filetype, 'cli')
    elseif strcmp(filetype, 'sp')
    end
    legend('Location', 'northeast');

    % Remove X Tick Labels except for the last subplot
    if k < nsubplots
        set(gca,'XTickLabel',[])
    end

    % Add stimulation marks
    ax = gca;
    line([stim_start_plot, stim_start_plot], [ax.YLim(1), ax.YLim(2)], ...
        'Color', 'r', 'LineStyle', '--');            % line for stimulation on
    line([stim_start_plot + stim_dur_plot, stim_start_plot + stim_dur_plot], [ax.YLim(1), ax.YLim(2)], ...
        'Color', 'r', 'LineStyle', '--');            % line for stimulation off

    % Add a stimulation line & title for the first subplot and an x-axis label for the last subplot
    if k == 1
        if strcmp(filetype, 'v')
            title(['Somatic voltage (mV) for ', strrep(filename, '_', '\_')]);
        elseif strcmp(filetype, 'cli')
            title(['Chloride concentration (mM) for ', strrep(filename, '_', '\_')]);
        elseif strcmp(filetype, 'sp')
            title(['Traces for ', strrep(filename, '_', '\_')]);
        end
    elseif k == nsubplots
        xlabel(timelabel);
    end
end

% Save figure
save_all_figtypes(h, figname, figtypes);
%close(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [timevec, timelabel, xlim1, xlim2, stim_start_plot, stim_dur_plot, spiketimevec] = set_time_units(data, tstart, tstop, stim_start, stim_dur, spikes)

% Change units of time axis from ms to s if the total time is > 10 seconds
if tstop > 10000
    timevec = data(:, 1)/1000;
    timelabel = 'Time (s)';
    xlim1 = tstart/1000;
    xlim2 = tstop/1000;
    stim_start_plot = stim_start/1000;
    stim_dur_plot = stim_dur/1000;
    if nargin >= 6
        spiketimevec = spikes(:, 2)/1000;
    else
        spiketimevec = [];
    end
else
    timevec = data(:, 1);
    timelabel = 'Time (ms)';
    xlim1 = tstart;
    xlim2 = tstop;
    stim_start_plot = stim_start;
    stim_dur_plot = stim_dur;
    if nargin >= 6
        spiketimevec = spikes(:, 2);
    else
        spiketimevec = [];
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}
