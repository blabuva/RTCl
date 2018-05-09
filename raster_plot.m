function [spikes, numactive, latency, actvel] = raster_plot(infolder, varargin)
%% Shows a spike raster plot and compute numactive, latency, actvel for each set of neurons (each .spi file in the infolder)
% USAGE: [spikes, numactive, latency, actvel] = raster_plot(infolder, varargin)
% Arguments:
%       infolder    - the name of the directory containing the .syn files, e.g. '20170317T1127_Ggaba_0.01'
%                   must be a directory
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
%                   - 'PlotSpikes': whether to plot raster plots
%                   must be logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PlotTuning': whether to plot tuning curves
%                   must be logical 1 (true) or 0 (false)
%                   default == true
%
% Requires:
%        cd/tuning_curves.m
%        cd/tuning_maps.m
%        infolder/*.spi
%        infolder/*loopedparams.mat
%        infolder/['sim_params_', pstring, '.csv'] 
%           for all the possible parameter strings
%        /home/Matlab/Downloaded_Functions/dirr.m
%        /home/Matlab/Downloaded_Functions/subaxis.m
%        /home/Matlab/Adams_Functions/isfigtype.m
%        /home/Matlab/Adams_Functions/save_all_figtypes.m
%        /home/Matlab/Adams_Functions/plot_tuning_curve.m 
%                                                   (through tuning_curves.m)
%        /home/Matlab/Adams_Functions/plot_tuning_map.m 
%                                                   (through tuning_maps.m)
%        /home/Matlab/Adams_Functions/find_ind_str_in_cell.m
%        /home/Matlab/Adams_Functions/get_loopedparams.m
% Used by:
%        cd/neuronlaunch.m
%        cd/tuning_curves.m
%        cd/tuning_maps.m
%
% 2017-02-08 Modified from rate.m of MCM 201505 for Peter's project
% 2017-03-05 Made tstop and ncells arguments
% 2017-03-13 Added line for stimulation
% 2017-03-13 Mark stimulation with dotted lines instead
% 2017-03-13 Plot spikes with a line rather than with a dot
% 2017-03-14 Added tstart
% 2017-03-16 Time unit labels are now flexible
% 2017-03-21 Fixed stim_start and stim_dur when time units are changed and added stim_freq to stim label
% 2017-03-23 Added actcellID as an argument and numactive, latency, actvel as outputs
% 2017-03-27 Filenames are now modified in neuronlaunch.m
% 2017-03-30 Files are now read in original order
% 2017-03-31 Now reads simulation parameters from the csv file
% 2017-03-31 Fixed ylimits of tuning curves
% 2017-04-02 Changed colormap from lines to jet
% 2017-04-04 Extract pnames, plabels, pmin, pmax, pinc, pislog, ncells from a saved .mat file
% 2017-04-06 Added clear all and close all inside parfor
% 2017-04-10 Added renewparpool
% 2017-04-10 Made renewparpool, maxnumworkers arguments
% 2017-04-14 Moved tuning curves part to tuning_curves.m
% 2017-04-18 Added 'FigType' as a parameter-value pair argument
% 2017-05-01 Changed all optional arguments into parameter-value pairs
% 2017-05-01 Added stimcellIDs; separated spiketimes to spiketimes_stim and spiketimes_nonstim
% 2017-05-04 Changed definition of ntrials
% 2017-05-04 Expanded current parameter string to accept multiple parameters changed at a time
% 2017-05-04 Added loopmode
% 2017-05-04 latency_cell_to_plot not given anymore because it is dependent on actmode
% 2017-05-09 Changed PlotCurves -> PlotTuning
% 2017-05-09 Changed 'FigType' to 'FigTypes', which now accepts a cell array of figure types
% 2017-05-09 Replaced saveas() with save_all_figtypes()

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
addpath(fullfile(functionsdirectory, '/Downloaded_Functions/'));    % for dirr.m & subaxis.m
addpath(fullfile(functionsdirectory, '/Adams_Functions/'));        % for isfigtype.m, find_ind_str_in_cell.m 
                                    % & get_loopedparams.m

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error('An infolder is required, type ''help show_RTnet'' for usage');
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
addParameter(iP, 'PlotSpikes', true, ...        % whether to plot raster plots
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotTuning', true, ...        % whether to plot tuning curves
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, infolder, varargin{:});
[~, figtypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);
outfolder = iP.Results.OutFolder;
maxnumworkers = iP.Results.MaxNumWorkers;
renewparpool = iP.Results.RenewParpool;
plotspikes = iP.Results.PlotSpikes;
plottuning = iP.Results.PlotTuning;

% Change default arguments if necessary
if strcmp(outfolder, '@infolder')
    outfolder = infolder;
end
if ~plotspikes
    renewparpool = false;        % only plotting spikes will take up significant amount of memory
end

%% Find all .spi files
files = dirr(infolder, '.spi');
nfiles = length(files);

%% Get loop parameters
[nump, pnames, plabels, pislog, pvalues, nperp, pchnames, pchvalues, ncells, actmode, loopmode] = get_loopedparams(infolder);

%% Check number of .spi files
ntrials = numel(pchnames);
if nfiles ~= ntrials
    error('Number of .spi files and number of trials don''t match!');
end

%% Create raster plot
spikes = cell(ntrials, 1);        % raw spike data
numactive = zeros(ntrials, 1);        % number of cells activated
latency = zeros(ntrials, ncells);    % latency to activation (seconds) for each cell
actvel = zeros(ntrials, 1);        % activation velocity (cells/second)
ct = 0;                    % counts number of trials completed
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
while ct < ntrials            % while not trials are completed yet
    first = ct + 1;                % first trial in this batch
    if renewparpool && ct + numworkers <= ntrials    % if memory is to be released
        last = ct + numworkers;            % limit the batch to numworkers
    else
        last = ntrials;
    end
    if renewparpool
        poolobj = parpool('local', numworkers);    % recreate a parallel pool object using fewer workers to prevent running out of memory
    end
    parfor i = first:last
    %for i = first:last
        % Construct current parameter string
        pstring = '';            % initialize for parfor
        if iscell(pchvalues)
            pstring = construct_suffix('NameValuePairs', {pchnames{i}, pchvalues{i}});
        elseif isnumeric(pchvalues)
            pstring = construct_suffix('NameValuePairs', {pchnames{i}, pchvalues(i)});
        end
        pspifile = [pstring, '.spi'];

        % Load spike data
        jnow = 0;            % initialize for parfor
        for j = 1:nfiles
            if ~isempty(strfind(files(j).name, pspifile))
                jnow = j;    % index in files of current file
                spikes{i} = load(fullfile(infolder, files(jnow).name));
            end
        end

        % Don't plot if file not found or no spikes
        stimcellIDs = 0;        % initialize for parfor    
        if isempty(spikes{i})
            fprintf('The file %s does not exist or has no spikes!\n', pspifile);
            numactive(i) = 0;            % number of cells activated is zero
            latency(i, :) = Inf*ones(1, ncells);    % latency to activation (seconds) is infinite for every cell
            actvel(i) = 0;                % activation velocity (cells/second) is zero
        else
            % Extract parameters from sim_params file
            fid = fopen(fullfile(infolder, ['sim_params_', pstring, '.csv']));
            simfilecontent = textscan(fid, '%s %f %s', 'Delimiter', ',');
            paramnames = simfilecontent{1};
            params_val = simfilecontent{2};
            tstart = params_val(find_ind_str_in_cell('tstart', paramnames, 'SearchMode', 'exact'));
            tstop = params_val(find_ind_str_in_cell('tstop', paramnames, 'SearchMode', 'exact'));
            RERErad = params_val(find_ind_str_in_cell('RERErad', paramnames, 'SearchMode', 'exact'));
            actcellID = params_val(find_ind_str_in_cell('actcellID', paramnames, 'SearchMode', 'exact'));
            stim_start = params_val(find_ind_str_in_cell('stim_start', paramnames, 'SearchMode', 'exact'));
            stim_dur = params_val(find_ind_str_in_cell('stim_dur', paramnames, 'SearchMode', 'exact'));
            stim_freq = params_val(find_ind_str_in_cell('stim_freq', paramnames, 'SearchMode', 'exact'));
            fclose(fid);

            % Extract info
            spikecelln = spikes{i}(:, 1);
            spiketimes = spikes{i}(:, 2);
            xlim1 = tstart;
            xlim2 = tstop;
            stim_start_plot = stim_start;
            stim_dur_plot = stim_dur;
            timelabel = 'Spike time (ms)';
    
            % Find number of cells activated
            numactive(i) = length(unique(spikecelln));    % number of unique cells with spike times

            % Find latency to activation (seconds) for each cell
            latency_this = zeros(1, ncells);        % for parfor
            for j = 1:ncells        % for each cellID == j-1
                ind = find(spikecelln == j-1);
                if isempty(ind)    
                    latency_this(j) = Inf;        % make latency infinite if cell not activated
                else
                    latency_this(j) = (min(spiketimes(ind)) - stim_start) / 1000;
                end
            end
            latency(i, :) = latency_this;            % for parfor

            % Find activation velocity (cells/second) for single-neuron or 3-neuron activation
            if (actmode == 1 || actmode == 3 || actmode == 4)
                celln_far = min(spikecelln);            % the cellID of the farthest activated cell to the left
                % fprintf('celln_far for file #%d = %d\n', i, celln_far);
                if celln_far == actcellID || (actmode == 3 && celln_far >= actcellID-RERErad-1 )
                    actvel(i) = 0;                % set activation velocity to zero if no spread to the left
                else
                    % Set up vectors for regression
                    if actmode ~= 3
                        celln_reg = celln_far:(actcellID-1);        % the cellIDs to regress
                    else
                        celln_reg = celln_far:(actcellID-RERErad-2);        % the cellIDs to regress
                    end
                    cellind_reg = celln_reg + 1;            % the cell indices to regress
                    lat_reg = latency_this(cellind_reg);         % the corresponding activation latencies

                    % Remove infinite entries
                    celln_reg2 = celln_reg(~isinf(lat_reg));
                    lat_reg2 = lat_reg(~isinf(lat_reg));

                    % Normalize vectors and find activation velocity with linear regression
                    celln_reg3 = celln_reg2(1) - celln_reg2;    % the normalized cellIDs to regress
                    lat_reg3 = lat_reg2 - lat_reg2(1);        % the normalized activation latencies (seconds)
                    actvel(i) = regress(celln_reg3', lat_reg3');    % the regressed coefficient will be 
                                            %    the activation velocity (cells/second)

                    %{
                    % Find activation velocity with linear approximation at the ends
                    cellind_far = celln_far + 1;            % the index for the farthest activated cell to the left
                    celln_near = max(actcellID-1, celln_far);    % the cellID of the nearest activated cell to the left
                    cellind_near = celln_near + 1;            % the index for the nearest activated cell to the left
                    actvel(i) = (celln_near - celln_far) / ...
                            (latency_this(cellind_far) - latency_this(cellind_near));
                                            % The approximate activation velocity
                    %}
                end
            end    

            % Change units of time axis from ms to s if the total time is > 10 seconds
            if tstop > 10000
                spiketimes = spiketimes/1000;
                timelabel = 'Spike time (s)';
                xlim1 = xlim1/1000;
                xlim2 = xlim2/1000;
                stim_start_plot = stim_start/1000;
                stim_dur_plot = stim_dur/1000;
            end

            % Find maximum and minimum time to plot
            xmin = min([xlim1, min(spiketimes)]);
            xmax = max([xlim2, max(spiketimes)]);

            % Find the IDs of cells that are stimulated or artificially activated
            switch actmode
            case {1, 4}
                stimcellIDs = actcellID;
            case 2
                ct = 0;
                for c = 0:ncells-1
                    if mod(c - actcellID, RERErad + 1) == 0
                        stimcellIDs(ct) = actcellID;
                        ct = ct + 1;
                    end
                end
            case 3
                stimcellIDs = [actcellID, actcellID + RERErad + 1, actcellID - RERErad - 1];
            end

            % Separate the spike times and cell numbers into two vectors
            Lia1 = ismember(spikecelln, stimcellIDs);    % whether the index in spikecelln belongs to stimcellIDs
            Lia2 = (spiketimes <= stim_start_plot + stim_dur_plot) ...
                & (spiketimes >= stim_start_plot);    % whether the spike time is within the stimulation period
            spiketimes_stim = spiketimes(Lia1 & Lia2);
            spiketimes_nonstim = spiketimes(~(Lia1 & Lia2));
            spikecelln_stim = spikecelln(Lia1 & Lia2);
            spikecelln_nonstim = spikecelln(~(Lia1 & Lia2));

            % Plot raster plot
            if plotspikes
                % Create plot
            %    h = figure(floor(rand()*10^4+1));
                h = figure(30000);
                clf(h);
                hold on;

                % Create raster plot
                line([spiketimes_stim, spiketimes_stim]', ...
                    [spikecelln_stim - 0.5, spikecelln_stim + 0.5]', ...
                    'Color', 'c', 'LineWidth', 0.5);        % plot spikes from stimulated cells
                line([spiketimes_nonstim, spiketimes_nonstim]', ...
                    [spikecelln_nonstim - 0.5, spikecelln_nonstim + 0.5]', ...
                    'Color', 'b', 'LineWidth', 0.5);        % plot spikes from nonstimulated cells
                line([stim_start_plot, stim_start_plot], [-1, ncells], ...
                    'Color', 'r', 'LineStyle', '--');            % line for stimulation on
                text(stim_start_plot + 0.5, ncells*0.975, ['Stim ON: ', num2str(stim_freq), ' Hz'], ...
                    'Color', 'r');                        % text for stimulation on
                text(xlim1 + 0.5, ncells*0.925, ['Number of cells activated: ', num2str(numactive(i))], ...
                    'Color', 'k');                        % text for number of cells activated
                text(xlim1 + 0.5, ncells*0.875, ...
                    ['Latency to activation for Cell#', num2str(actcellID-1), ': ', ...
                    num2str(latency_this(actcellID-1+1)), ' s'], ...
                    'Color', 'k');                        % text for latency to activation (seconds)
                text(xlim1 + 0.5, ncells*0.825, ...
                    ['Activation velocity: ', num2str(actvel(i)), ' cells/sec'], ...
                    'Color', 'k');                        % text for activation velocity (cells/second)
                line([stim_start_plot + stim_dur_plot, stim_start_plot + stim_dur_plot], ...
                    [-1,  ncells], 'Color', 'r', 'LineStyle', '--');    % line for stimulation off
                text(stim_start_plot + stim_dur_plot + 0.5, ncells*0.975, 'Stim OFF', ...
                    'Color', 'r');                        % text for stimulation off
                xlim([xmin, xmax]);                        % time range to plot
                ylim([-1, ncells]);                        % cell ID runs from 0 to ncells-1
                xlabel(timelabel);
                ylabel('Neuron number');
                title(strrep(files(jnow).name, '_', '\_'));

                % Save figure
                figname = fullfile(outfolder, strrep(files(jnow).name, '.spi', '_raster_plot.png'));
                save_all_figtypes(h, figname, figtypes);
                % close(h);
            end    
        end
%        close all;
    end

    if renewparpool
        delete(poolobj);            % delete the parallel pool object to release memory
    end
    ct = last;                % update number of trials completed
end
if renewparpool
    poolobj = parpool('local', oldnumworkers);    % recreate a parallel pool object using the previous number of workers
end

%% Plot tuning curves
if plottuning
    switch loopmode
    case 'cross'
        latency_cells_to_plot = 0:ncells/2-1;        % cells to plot for the latency tuning curve
        tuning_curves (infolder, outfolder, numactive, latency, actvel, ...
                nump, pnames, plabels, pislog, pvalues, nperp, ...
                ncells, actmode, loopmode, 'CellsToPlot', latency_cells_to_plot, ...
                'FigTypes', figtypes);
    case 'grid'
        tuning_maps (infolder, outfolder, numactive, latency, actvel, ...
                nump, pnames, plabels, pislog, pvalues, nperp, ...
                ncells, actmode, loopmode, 'FigTypes', figtypes);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}


