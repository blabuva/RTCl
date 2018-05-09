%% neuronlaunch.m
%% Launches NEURON with simulation commands and plot output figures
%
% Requires:
% Dependent files in cd/:
%       cd/raster_plot.m            (through neuronlaunch.m, tuning_curves.m,
%                                       tuning_maps.m)
%       cd/show_RTnet.m             (through neuronlaunch.m)
%       cd/single_neuron.m          (through neuronlaunch.m)
%       cd/tuning_curves.m          (through raster_plot.m)
%       cd/tuning_maps.m            (through raster_plot.m)
%       cd/run.hoc                  (through neuronlaunch.m)
%       cd/run_small.hoc            (through neuronlaunch.m)
%       cd/net.hoc                  (through run.hoc, run_small.hoc)
%       cd/RE.tem                   (through net.hoc)
%       cd/cadecay.mod              (through RE.tem)
%       cd/cldecay.mod              (through RE.tem)
%       cd/cldecay1.mod             (through RE.tem)
%       cd/cldecay2.mod             (through RE.tem)
%       cd/gabaA_Cl.mod             (through net.hoc)
%       cd/HH2.mod                  (through RE.tem)
%       cd/ipulse2.mod              (through net.hoc)
%       cd/ITs.mod                  (through RE.tem)
%       cd/IKCa.mod                 (through RE.tem)
%
% Dependent files in functionsdirectory/Adams_Functions/:
%       all_ordered_pairs.m         (through make_loopedparams.m, 
%                                       get_loopedparams.m)
%       change_params.m             (through neuronlaunch.m)
%       check_and_collapse_identical_contents.m
%                                   (through tuning_curves.m)
%       construct_fullfilename.m    (through neuronlaunch.m)
%       construct_suffix.m          (through construct_fullfilename.m)
%       find_ind_str_in_cell.m      (through neuronlaunch.m,
%                                       raster_plot.m, single_neuron.m,
%                                       tuning_curves.m, tuning_maps.m,
%                                       validate_string.m)
%       get_loopedparams.m          (through raster_plot.m, 
%                                       tuning_curves.m, tuning_maps.m)
%       intersectm.m                (through find_ind_str_in_cell.m)
%       isfigtype.m                 (through raster_plot.m, single_neuron.m,
%                                       tuning_curves.m, tuning_maps.m,
%                                       plot_tuning_curve.m, plot_tuning_map.m,
%                                       save_all_figtypes.m)
%       make_loopedparams.m         (through neuronlaunch.m)
%       plot_tuning_curve.m         (through raster_plot.m, tuning_curves.m)
%       plot_tuning_map.m           (through raster_plot.m, tuning_maps.m)
%       save_all_figtypes.m         (through raster_plot.m, single_neuron.m,
%                                       tuning_curves.m, tuning_maps.m,
%                                       plot_tuning_curve.m, plot_tuning_map.m)
%       update_params.m             (through change_params.m)
%       validate_string.m           (through isfigtype.m)
%       vec2cell.m                  (through tuning_curves.m)
%       vec2array.m                 (through tuning_maps.m)
%
% Dependent files in functionsdirectory/Downloaded_Functions/:
%       dirr.m                      (through show_RTnet.m, raster_plot.m,
%                                       single_neuron.m, tuning_curves.m)
%       subaxis.m                   (through raster_plot.m, single_neuron.m)
%       parseArgs.m                 (through subaxis.m)
%
% 2017-02-07 Modified from neuronlaunch.m of MCM 201505 for Peter's project
% 2017-02-24 Moved parameters and simulation commands 
%               from run.hoc to neuronlaunch.m
% 2017-02-24 Changed output file names
% 2017-02-26 Added chloride concentration traces
% 2017-03-04 Added spcellID, sp_data, etc.
% 2017-03-05 Let the variables stabilize for 3000 ms before current injection 
%               (chenged cp_start & tstop)
% 2017-03-07 Now records two special neurons at once
% 2017-03-14 Change current pulse stimulation to Ipulse2 from IClamp1
% 2017-03-14 Renamed tstart->timer; tstart is now time to start plotting
% 2017-03-14 Removed gaba_grel as a parameter and added useca, cldnum, tauKCC2
% 2017-03-14 Renamed RErest->REepas; make that an argument for init() in RE
% 2017-03-15 Added REuseca to make T channels, KCa channels and 
%               calcium decay optional
% 2017-03-15 Added drive_ch, drive_ex, cli1
% 2017-03-16 Added simmode
% 2017-03-21 Changed stimulation time to 10 seconds and 
%               relaxation time to 10 seconds
% 2017-03-21 Added REnsegs
% 2017-03-22 Added REconsyn
% 2017-03-22 Changed order of proplabels
% 2017-03-27 Added parameters to loop through and print and save everything 
%               separately for each iteration
% 2017-03-30 Made cp_amp dependent on RE_diam^2
% 2017-03-31 Added trialn to save
% 2017-03-31 Changed units of REtauKCC2 to seconds
% 2017-04-04 Moved IDlabels and proplabels to single_neuron.m
% 2017-04-05 Added the actmode (now 2) for REmultcp()
% 2017-04-06 Do not save spikes and data in matfile
% 2017-04-10 Added the actmode (now 3) for REthreecp()
% 2017-04-10 Added renewparpoolflag
% 2017-04-10 Changed cp_num from floor(stim_dur/cp_per) to ceil(stim_dur/cp_per)
% 2017-04-17 Changed ntrperp to nperp
% 2017-04-17 Now generates and saves pvalues and nperp as in get_loopparams.m
% 2017-05-02 Added onlargememflag
% 2017-05-03 Added singletrialnum to run just one trial
% 2017-05-03 Now uses pchnames and pchvalues directly to change parameters
% 2017-05-03 Added loopmode and implemented the 'grid' case
% 2017-05-04 Changed argument scheme for construct_fullfilename()
% 2017-05-08 Changed onlargememflag to onRivannaflag
% 2017-05-08 Revert back to onlargememflag
% 2017-05-09 Changed plotcurves to plottuning
% 2018-05-09 Found all dependent files and listed them
% 2018-05-09 Now compiles .mod files at the beginning of each run

%% Clear workspace and close figures
clear all force hidden
close all force hidden

%% Compile custom .mod files
unix('nrnivmodl');

%% This is for MATLAB R2017a and beyond, for compatibility with the code
set(0, 'defaultLegendAutoUpdate', 'off');

%% Experiment Name
experimentname = 'RTCl';

%% Flags
debugflag = 0;              % very short simulation
singletrialnum = 0;         % run only one trial with this trial number
onlargememflag = 0;         % run on large memory nodes
saveplotmode = 'curves';    % what to save and plot
loopmode = 'cross';         % how to loop through parameters: 'cross' or 'grid'
                            %   'cross' - Loop through each parameter 
                            %               while fixing others
                            %   'grid'  - Loop through all possible 
                            %               combinations of parameters

%% For parpool
if onlargememflag           % using Rivanna
    renewparpoolflag_NEURON = 0;    % whether to renew parpool every batch 
                                    %   to release memory for NEURON simulations
    maxnumworkers_NEURON = 20;      % maximum number of workers 
                                    %   for running NEURON 
    renewparpoolflag_plots = 0;     % whether to renew parpool every batch 
                                    %   to release memory for MATLAB plotting
    maxnumworkers_plots = 20;       % maximum number of workers 
                                    %   for MATLAB plotting
else
    switch saveplotmode
    case 'all'              % saving and plotting everything
        renewparpoolflag_NEURON = 1;
        maxnumworkers_NEURON = 8;
        renewparpoolflag_plots = 1;
        maxnumworkers_plots = 8;
    case 'curves'           % saving spikes and plotting curves/maps only
        renewparpoolflag_NEURON = 0;
        maxnumworkers_NEURON = 20;
        renewparpoolflag_plots = 0;
        maxnumworkers_plots = 20;
    case 'spikes'           % saving spikes and plotting raster plots 
                            %   and curves/maps only
        renewparpoolflag_NEURON = 0
        maxnumworkers_NEURON = 20; 
        renewparpoolflag_plots = 1;
        maxnumworkers_plots = 10; 
    case 'spikes&special'   % saving spikes and special neuron traces only
        renewparpoolflag_NEURON = 0;
        maxnumworkers_NEURON = 12;
        renewparpoolflag_plots = 1;
        maxnumworkers_plots = 10;
    end
end

switch saveplotmode
case 'all'
    %% Save flags
    savenetwork = 1;        % whether to save network topology
    savespikes = 1;         % whether to save spike data
    savesomavoltage = 1;    % whether to save all voltage data
    savesomacli = 1;        % whether to save all chloride concentration data
    savespecial = 1;        % whether to save special neuron data

    %% Plot flags
    plotnetwork = 1;        % whether to plot network topology
    plotspikes = 1;         % whether to plot spike data
    plottuning = 1;         % whether to plot tuning curves
    plotsingleneurondata = 1;   % whether to plot single neuron data
case 'curves'
    %% Save flags
    savenetwork = 0;        % whether to save network topology
    savespikes = 1;         % whether to save spike data
    savesomavoltage = 0;    % whether to save all voltage data
    savesomacli = 0;        % whether to save all chloride concentration data
    savespecial = 0;        % whether to save special neuron data

    %% Plot flags
    plotnetwork = 0;        % whether to plot network topology
    plotspikes = 0;         % whether to plot spike data
    plottuning = 1;         % whether to plot tuning curves
    plotsingleneurondata = 0;   % whether to plot single neuron data
case 'spikes'
    %% Save flags
    savenetwork = 0;        % whether to save network topology
    savespikes = 1;         % whether to save spike data
    savesomavoltage = 0;    % whether to save all voltage data
    savesomacli = 0;        % whether to save all chloride concentration data
    savespecial = 0;        % whether to save special neuron data

    %% Plot flags
    plotnetwork = 0;        % whether to plot network topology
    plotspikes = 1;         % whether to plot spike data
    plottuning = 1;         % whether to plot tuning curves
    plotsingleneurondata = 0;   % whether to plot single neuron data
case 'spikes&special'
    %% Save flags
    savenetwork = 0;        % whether to save network topology
    savespikes = 1;         % whether to save spike data
    savesomavoltage = 0;    % whether to save all voltage data
    savesomacli = 0;        % whether to save all chloride concentration data
    savespecial = 1;        % whether to save special neuron data

    %% Plot flags
    plotnetwork = 0;        % whether to plot network topology
    plotspikes = 1;         % whether to plot spike data
    plottuning = 1;         % whether to plot tuning curves
    plotsingleneurondata = 1;   % whether to plot single neuron data
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters to loop through
%{
pnames  = {'REdiam'};       % names of parameters to loop through
plabels = {'REdiam (um)'};  % labels of parameters to loop through
pmin    = [2];              % minimum values of parameters to loop through
pmax    = [20];             % maximum values of parameters to loop through
pinc    = [2];              % increments of parameters to loop through
pislog  = [0];              % whether increments of parameters is in log
%}
%{
pnames  = {'REGgaba'};      % names of parameters to loop through
plabels = {'REGgaba (uS)'}; % labels of parameters to loop through
pmin    = [0.0025];         % minimum values of parameters to loop through
pmax    = [0.045];          % maximum values of parameters to loop through
pinc    = [0.0025];         % increments of parameters to loop through
pislog  = [0];              % whether increments of parameters is in log
%}
%{
pnames  = {'stim_freq'};    % names of parameters to loop through
plabels = {'Stimulation Frequency (Hz)'};   
                            % labels of parameters to loop through
pmin    = [1];              % minimum values of parameters to loop through
pmax    = [128];            % maximum values of parameters to loop through
pinc    = [2^(1/2)];        % increments of parameters to loop through
pislog  = [1];              % whether increments of parameters is in log
%}
%{
pnames  = {'REtauKCC2'};    % names of parameters to loop through
plabels = {'Time constant of KCC2 (s)'};    
                            % labels of parameters to loop through
pmin    = [4];              % minimum values of parameters to loop through
pmax    = [64];             % maximum values of parameters to loop through
pinc    = [2^(1/4)];        % increments of parameters to loop through
pislog  = [1];              % whether increments of parameters is in log
%}
%{
pnames  = {'REdiam', 'REGgaba'};
                            % names of parameters to loop through
plabels = {'REdiam (um)', 'REGgaba (uS)'};
                            % labels of parameters to loop through
pmin    = [8, 0.1];         % minimum values of parameters to loop through
pmax    = [12, 0.5];        % maximum values of parameters to loop through
pinc    = [0.5, 0.05];      % increments of parameters to loop through
pislog  = [0, 0];           % whether increments of parameters is in log
%}

%% Final parameters to loop through

%{
% neuronlaunch99~101
pnames  = {'stim_freq', 'REtauKCC2'};
                            % names of parameters to loop through
plabels = {'Stimulation Frequency (Hz)', 'Time constant of KCC2 (s)'};
                            % labels of parameters to loop through
pmin    = [2, 4];           % minimum values of parameters to loop through
pmax    = [58, 64];         % maximum values of parameters to loop through
pinc    = [4, 2^(1/4)];     % increments of parameters to loop through
pislog  = [0, 1];           % whether increments of parameters is in log
%}
%{
% neuronlaunch102.slurm
pnames  = {'stim_freq', 'REtauKCC2'};
                            % names of parameters to loop through
plabels = {'Stimulation Frequency (Hz)', 'Time constant of KCC2 (s)'};
                            % labels of parameters to loop through
pmin    = [32, 32];         % minimum values of parameters to loop through
pmax    = [52, 64];         % maximum values of parameters to loop through
pinc    = [0.2, 2^(1/64)];  % increments of parameters to loop through
pislog  = [0, 1];           % whether increments of parameters is in log
stim_freq = 14;             % stimulation frequency (Hz), 
                            %   must be less than 1000/cp_dur
REtauKCC2 = 8;              % Cl- removal time constant (s) in RE cells
%}
%{
% neuronlaunch103.slurm
pnames  = {'stim_freq', 'REtauKCC2'};
                            % names of parameters to loop through
plabels = {'Stimulation Frequency (Hz)', 'Time constant of KCC2 (s)'};
                            % labels of parameters to loop through
pmin    = [20, 16];         % minimum values of parameters to loop through
pmax    = [35, 64];         % maximum values of parameters to loop through
pinc    = [0.2, 2^(1/64)];  % increments of parameters to loop through
pislog  = [0, 1];           % whether increments of parameters is in log
stim_freq = 20;             % stimulation frequency (Hz)
                            %   must be less than 1000/cp_dur
REtauKCC2 = 16;             % Cl- removal time constant (s) in RE cells
%}
%{
% neuronlaunch104.slurm
pnames  = {'stim_freq', 'REtauKCC2'};
                            % names of parameters to loop through
plabels = {'Stimulation Frequency (Hz)', 'Time constant of KCC2 (s)'};
                            % labels of parameters to loop through
pmin    = [14, 8];          % minimum values of parameters to loop through
pmax    = [28, 16*2^(1/2)]; % maximum values of parameters to loop through
pinc    = [0.2, 2^(1/64)];  % increments of parameters to loop through
pislog  = [0, 1];           % whether increments of parameters is in log
stim_freq = 30;             % stimulation frequency (Hz)
                            %   must be less than 1000/cp_dur
REtauKCC2 = 32;             % Cl- removal time constant (s) in RE cells
%}
%{
% neuronlaunch105~108
pnames  = {'stim_freq'};    % names of parameters to loop through
plabels = {'Stimulation Frequency (Hz)'};   
pmin    = [1];              % minimum values of parameters to loop through
pmax    = [38];             % maximum values of parameters to loop through
pinc    = [1];              % increments of parameters to loop through
pislog  = [0];              % whether increments of parameters is in log
%}
%{
% neuronlaunch109~111
pnames  = {'stim_freq'};    % names of parameters to loop through
plabels = {'Stimulation Frequency (Hz)'};   
pmin    = [10];             % minimum values of parameters to loop through
pmax    = [30];             % maximum values of parameters to loop through
pinc    = [0.1];            % increments of parameters to loop through
pislog  = [0];              % whether increments of parameters is in log
%}
%{
% neuronlaunch112
pnames  = {'stim_freq'};    % names of parameters to loop through
plabels = {'Stimulation Frequency (Hz)'};   
pmin    = [0.1];            % minimum values of parameters to loop through
pmax    = [20];             % maximum values of parameters to loop through
pinc    = [0.1];            % increments of parameters to loop through
pislog  = [0];              % whether increments of parameters is in log
%}

% neuronlaunch109~111
pnames  = {'stim_freq'};    % names of parameters to loop through
plabels = {'Stimulation Frequency (Hz)'};   
pmin    = [10];             % minimum values of parameters to loop through
pmax    = [30];             % maximum values of parameters to loop through
pinc    = [5];              % increments of parameters to loop through
pislog  = [0];              % whether increments of parameters is in log

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Global parameters to be defined at the start of NEURON, 
%   to be consistent with run.hoc or run.hoc
ncells = 10; %100;       % SET IN run.hoc: total number of cells in the network
celsius = 31;       % SET IN run.hoc: temperature of experiment (celsius), 
                    %                   Peter's value 
                    %                   (Sohal & Huguenard 2003 used 34 degC)
nsp = 2;            % SET IN run.hoc: number of special neurons
RERErad = 4;        % SET IN run.hoc: radius of intra-RE connections,
                    %                   Sohal & Huguenard 2003

%% Network parameters
sp_thr = 0;         % action potential threshold (mV)
syn_del = 1;        % synaptic delay (ms)
syn_w = 1;          % synaptic weight (fraction of channels activated)
                    %   for simplicity assume channels are always activated 
                    %   and that channels have linearly additive effects

%% RE cell parameters
REnsegs = 1;        % number of segments in an RE cell (1, 3 or 9)
                    %   if REnsegs >= 3, the GABAA synapses will be 
                    %   distributed on each side
REuseca = 2;        % mode to use ca_ion (0 or 1 or 2) in RE cells
                    %   0 - gTs & gKCa set to zero for all cells
                    %   1 - Calcium mechanisms inserted and normal in all cells
                    %   2 - gTs & gKCa set to zero for stimulated cells
REcldnum = 2;       % which cld mechanism to use (0 or 1 or 2) in RE cells
REconsyn = 0;       % whether to concentrate synapses (0 or 1) in RE cells
                    %   for REnsegs = 1 or 3, 
                    %   we have soma, soma_flank[0], soma_flank[1]
REtauKCC2 = 32;     % Cl- removal time constant (s) in RE cells
REepas = -70;       % leak reversal potential (mV) of RE cells, Peter's value
                    %   Sohal & Huguenard 2003 used -77 mV
                    %   Jedlicka et al 2011 used -60 mV
REdiam = 10;        % diameter (um) of an RE cell, Peter's value
REgpasLB = 4.5e-5;  % lower bound for passive leak conductance (S/cm^2) 
                    %   in RE cells, Sohal & Huguenard 2003
REgpasUB = 5.5e-5;  % upper bound for passive leak conductance (S/cm^2) 
                    %   in RE cells, Sohal & Huguenard 2003
                    %   Jedlicka et al 2011 used 2e-4 S/cm^2
REGgaba = 0.0025;   % peak conductance (uS) for each event 
                    %   at the GABA-A receptors 
                    %   (Sohal & Huguenard 2003 varied between 5~12.5 nS)

%% Initial ion concentrations
cai0 = 2.4e-4;      % initial intracellular [Ca++] (mM), Destexhe et al
cao0 = 2;           % initial extracellular [Ca++] (mM), Peter's value
% cli0 = 5;         % initial intracellular [Cl-] (mM), 
                    %   corresponding to eGABA = -71 mV
                    %            TODO: Ulrich
cli0 = 8;           % initial intracellular [Cl-] (mM), 
                    %   corresponding to eGABA = -61 mV
                    %            Jedlicka et al 2011 (agrees with Peter's data)
% cli0 = 11         % initial intracellular [Cl-] (mM), 
                    %   corresponding to eGABA = -54 mV
                    %            TODO: Peter's data with ChABC added
% cli0 = 17         % initial intracellular [Cl-] (mM), 
                    %   corresponding to eGABA = -45 mV
                    %            TODO: Sun
clo0 = 130.5;       % initial extracellular [Cl-] (mM), 
                    %   Peter's value (Jedlicka et al 2011 used 133.5 mM)

%% Activation mode
actmode = 3;        % 1 - Activate a single RE cell by injecting 
                    %       a train of current pulses
                    % 2 - Activate every (RERErad + 1)th RE cell 
                    %       by injecting trains of current pulses
                    % 3 - Activate 3 RE cells by injecting trains 
                    %       of current pulses
                    % 4 - Activate a single RE cell by changing 
                    %       the membrane potential instantaneously
                    % 5 - Activate RE cells with a Gaussian likelihood 
                    %       by changing the mp instantaneously

%% Activation parameters for 'cp' mode
actcellID = floor(ncells/2);    % ID # of central neuron to activate
stim_start = 3000;          % stimulation delay (ms)
stim_dur = 200000;          % stimulation duration (ms)
stim_freq = 14;             % stimulation frequency (Hz), 
                            %   must be less than 1000/cp_dur
cp_dur = 0.1;               % current pulse duration (ms)
cp_amp = 4*(REdiam/10)^2;   % current pulse amplitude (nA), 
                            %   must be proportional to square of diameter 

% The following must be consistent with update_params.m
cp_per = floor(1000/stim_freq);     % current pulse period (ms), 
                                    %   i.e. interval between pulse onsets
cp_num = ceil(stim_dur/cp_per);     % number of current pulses

%% Activation parameters for 'single' mode
actcellv = 0;               % voltage (mV) to set activated neuron to

%% Activation parameters for 'random' mode
actwidth = floor(ncells/2); % width of Gaussian distribution for 
                            %   randomly activating cells
actmaxp = 0.5;              % maximum likelihood of activation at center

%% Simulation parameters
simmode = 1;                % 1 - full simulation
                            % 2 - short simulation
                            % 3 - medium simulation
if simmode == 1
    tstop = 233000;         % total time of simulation (ms)
elseif simmode == 2
    tstop = 4000;           % total time of simulation (ms)
elseif simmode == 3
    tstop = 10000;          % total time of simulation (ms)
end
dt = 0.1;                   % time step of integration (ms)

%% Recording parameters
if simmode == 1
    tstart = 0;             % time to start plotting (ms)
elseif simmode == 2
    tstart = 2900;          % time to start plotting (ms)
elseif simmode == 3
    tstart = 2000;          % time to start plotting (ms)
end
sp1cellID = actcellID;      % ID # of 1st special neuron to record
sp2cellID = actcellID - 1;  % ID # of 2nd special neuron to record
if ncells == 100
    sp3cellID = actcellID - 10; 
                            % ID # of 3rd special neuron to record
elseif ncells == 10
    sp3cellID = actcellID - 5;  
                            % ID # of 3rd special neuron to record
end

%% Set ID #s of neurons to plot
act = actcellID;            % ID # of the activated neuron
act_left1 = actcellID - 1;  % ID # of the neuron one below the activated neuron
act_left2 = actcellID - 2;  % ID # of the neuron 2 below the activated neuron
far = actcellID - floor(ncells/10);
                            % ID # of a far away neuron
if ncells == 100
    far = actcellID - 10;
                            % ID # of a far away neuron
elseif ncells == 10
    far = actcellID - 5;
                            % ID # of a far away neuron
end

%% Arguments for plotting (not logged in sim_params)
properties = 1:1:12;        % property #s of special neuron to record 
                            %   to be plotted (maximum range: 1~12, 
                            %   must be consistent with net.hoc)
% properties = [1, 2, 3, 4, 8, 10, 11];
IDs = [act, act_left1, act_left2, far];
                            % ID #s for neurons whose voltage is to be plotted

%% Set output file names; must have only one '.' (not logged in sim_params)
simparamsF = 'sim_params.csv';  % file with simulation parameters
scmdsF = 'sim_commands.txt';    % file with simulation commands
soutF = 'sim_output.txt';       % file with simulation standard outputs
sREREsynF = 'RERE.syn';         % file with RE-RE synaptic connections
sREspikeF = 'RE.spi';           % file with RE spike train output
sREvF = 'RE.singv';             % file with RE single neuron voltage traces
sREcliF = 'RE.singcli';         % file with RE single neuron 
                                %   chloride concentration traces
sREsp1F = ['RE[', num2str(sp1cellID), '].singsp'];
                                % file with RE special neuron #1 other traces
sREsp2F = ['RE[', num2str(sp2cellID), '].singsp'];
                                % file with RE special neuron #2 other traces
sREsp3F = ['RE[', num2str(sp3cellID), '].singsp'];
                                % file with RE special neuron #3 other traces
sREleakF = 'REleak.csv';        % file with RE neuron leak properties
sREparamsF = 'REparams.csv';    % file with RE neuron parameters

%% For debug mode
if debugflag
%    dt = 0.1;
    tstart = 0;
    stim_start = 10;
    stim_dur = 100;
    tstop = 100;

    % Minimize number of points
    for p = 1:length(pnames)
        if pislog(p)
            pinc(p) = pmax(p)/pmin(p);
        else
            pinc(p) = pmax(p) - pmin(p);
        end
    end

    renewparpoolflag_NEURON = 0;
    maxnumworkers_NEURON = 12;
    renewparpoolflag_plots = 0;
    maxnumworkers_plots = 12;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set folders for reading and saving files
% Find home directory
if exist(fullfile(pwd, '/Adams_Functions/'), 'dir') == 7
    homedirectory = pwd;
elseif exist('/media/adamX/RTCl/', 'dir') == 7
    homedirectory = '/media/adamX/RTCl/';
elseif exist('/scratch/al4ng/RTCl/', 'dir') == 7
    homedirectory = '/scratch/al4ng/RTCl/';
else
    homedirectory = pwd;
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
addpath(fullfile(functionsdirectory, '/Adams_Functions/'));

% Make directory to save all data, 
%   use current date & time in the format: YYYYMMDDThhmm
ts1 = datestr(clock, 30);       % date & time in the format YYYYMMDDThhmmss
outfoldername = ts1(1:end-2);   % name for outfolder
outfolder = fullfile(homedirectory, outfoldername); % take off seconds
if exist(outfolder, 'dir') ~= 7
    mkdir(outfolder);
end

% Read data from same directory
infolder = outfolder;

%% Construct looped parameters
[pchnames, pchvalues, ntrials, nump, pvalues, nperp] = ...
    make_loopedparams (loopmode, pnames, plabels, pislog, pmin, pmax, pinc, ...
            'OutFolder', outfolder, 'FileLabel', outfoldername, ...
            'NCells', ncells, 'ActMode', actmode);

%% Create arrays for parameters
paramlabels = {
    '# of cells', 'temperature of experiment (celsius)', ...
    'number of special neurons', ...
    'radius of intra-RE connections', ...
    'action potential threshold (mV)', 'synaptic delay (ms)', ...
    'synaptic weight (fraction of channels activated)', ...
    'number of segments in an RE cell (must be odd)', ...
    'mode to use ca_ion (0 or 1 or 2) in RE cells', ...
    'which cld mechanism to use (0 or 1 or 2) in RE cells', ...
    'whether to concentrate synapses (0 or 1) in RE cells', ...
    'Cl- removal time constant (ms) in RE cells', ...
    'leak reversal potential (mV) of RE cells', ...
    'diameter (um) of an RE cell', ...
    'lower bound for passive leak conductance (S/cm^2) in RE cells', ...
    'upper bound for passive leak conductance (S/cm^2) in RE cells', ...
    'conductance (uS) of GABA-A synapses on RE cells', ...
    'initial intracellular [Ca++] (mM)', ...
    'initial extracellular [Ca++] (mM)', ...
    'initial intracellular [Cl-] (mM)', ...
    'initial extracellular [Cl-] (mM)', ...
    'activation mode', 'ID # of central neuron to activate', ...
    'stimulation delay (ms)', 'stimulation duration (ms)', ...
    'stimulation frequency (Hz)', ...
    'current pulse duration (ms)', 'current pulse amplitude (nA)', ...
    'current pulse period (ms)', 'number of current pulses', ...
    'voltage (mV) to set activated neuron to', ...
    'width of Gaussian distribution for randomly activating cells', ...
    'maximum likelihood of activation at center', ...
    'simulation mode', 'total time of simulation (ms)', ...
    'time step of integration (ms)', 'time to start plotting (ms)', ...
    'ID # of 1st special neuron to record', ...
    'ID # of 2nd special neuron to record', ...
    'ID # of 3rd special neuron to record', ...
    'whether to save network topology', 'whether to save spike data', ...
    'whether to save all voltage data', ...
    'whether to save all chloride concentration data', ...
    'whether to save special neuron data', ...
    'whether to plot network topology', ...
    'whether to plot spike data', ...
    'whether to plot single neuron data', ...
    'ID # of the activated neuron', ...
    'ID # of the neuron one below the activated neuron', ...
    'ID # of the neuron 2 below the activated neuron', ...
    'ID # of a far away neuron', ...
    'number of times to run simulation', 'current trial number'};
paramnames = { ...
    'ncells', 'celsius', 'nsp', 'RERErad', ...
    'sp_thr', 'syn_del', 'syn_w', ...
    'REnsegs', 'REuseca', 'REcldnum', ...
    'REconsyn', 'REtauKCC2', ...
    'REepas', 'REdiam', 'REgpasLB', 'REgpasUB', 'REGgaba', ...
    'cai0', 'cao0', 'cli0', 'clo0', ...
    'actmode', 'actcellID', ...
    'stim_start', 'stim_dur', 'stim_freq', ...
    'cp_dur', 'cp_amp', 'cp_per', 'cp_num', ...
    'actcellv', 'actwidth', 'actmaxp', ...
    'simmode', 'tstop', ...
    'dt', 'tstart', ...
    'sp1cellID', 'sp2cellID', ...
    'sp3cellID', ...
    'savenetwork', 'savespikes', 'savesomavoltage', ...
    'savesomacli', 'savespecial', ...
    'plotnetwork', 'plotspikes', 'plotsingleneurondata', ...
    'act', 'act_left1', 'act_left2', 'far', ...
    'ntrials', 'trialn'};
paramsinit = [ ...
    ncells, celsius, nsp, RERErad, ...
    sp_thr, syn_del, syn_w, ...
    REnsegs, REuseca, REcldnum, ...
    REconsyn, REtauKCC2, ...
    REepas, REdiam, REgpasLB, REgpasUB, REGgaba, ...
    cai0, cao0, cli0, clo0, ...
    actmode, actcellID, ...
    stim_start, stim_dur, stim_freq, ...
    cp_dur, cp_amp, cp_per, cp_num, ...
    actcellv, actwidth, actmaxp, ...
    simmode, tstop, ...
    dt, tstart, ...
    sp1cellID, sp2cellID, sp3cellID, ...
    savenetwork, savespikes, savesomavoltage, ...
    savesomacli, savespecial, ...
    plotnetwork, plotspikes, plotsingleneurondata, ...
    act, act_left1, act_left2, far, ...
    ntrials, singletrialnum];

%% Setup parameters and filenames for each trial
paramvals = cell(1, ntrials);            % stores parameters used for each trial 
simparamsF_full = cell(1, ntrials);
sREREsynF_full = cell(1, ntrials);
sREspikeF_full = cell(1, ntrials);
sREvF_full = cell(1, ntrials);
sREcliF_full = cell(1, ntrials);
sREsp1F_full = cell(1, ntrials);
sREsp2F_full = cell(1, ntrials);
sREsp3F_full = cell(1, ntrials);
sREleakF_full = cell(1, ntrials);
sREparamsF_full = cell(1, ntrials);
for k = 1:ntrials
    % Update trial count
    indtrialn = find_ind_str_in_cell('trialn', paramnames);
    paramvals{k}(indtrialn) = k;        % update trial number

    % Set parameters for this trial according to pchnames and pchvalues
    pchnames_k = pchnames{k};
    if iscell(pchvalues)
        pchvalues_k = pchvalues{k};
    elseif isnumeric(pchvalues)
        pchvalues_k = pchvalues(k);
    end
    [paramvals{k}] = change_params(pchnames_k, pchvalues_k, ...
                                    paramnames, paramsinit, ...
                                    'ExperimentName', experimentname);

    % Print parameters to a comma-separated-value file
    simparamsF_full{k} = construct_fullfilename(simparamsF, ...
                                'Directory', outfolder, ...
                                'NameValuePairs', {pchnames_k, pchvalues_k});
    fid = fopen(simparamsF_full{k}, 'w');
    for i = 1:length(paramvals{k})
        fprintf(fid, '%s, %g, %s\n', ...
            paramnames{i}, paramvals{k}(i), paramlabels{i});
    end
    fclose(fid);

    % Construct full file names
    sREREsynF_full{k}     = construct_fullfilename(sREREsynF, ...
                                'Directory', outfolder, ...
                                'NameValuePairs', {pchnames_k, pchvalues_k});
    sREspikeF_full{k}     = construct_fullfilename(sREspikeF, ...
                                'Directory', outfolder, ...
                                'NameValuePairs', {pchnames_k, pchvalues_k});
    sREvF_full{k}         = construct_fullfilename(sREvF, ...
                                'Directory', outfolder, ...
                                'NameValuePairs', {pchnames_k, pchvalues_k});
    sREcliF_full{k}     = construct_fullfilename(sREcliF, ...
                                'Directory', outfolder, ...
                                'NameValuePairs', {pchnames_k, pchvalues_k});
    sREsp1F_full{k}     = construct_fullfilename(sREsp1F, ...
                                'Directory', outfolder, ...
                                'NameValuePairs', {pchnames_k, pchvalues_k});
    sREsp2F_full{k}     = construct_fullfilename(sREsp2F, ...
                                'Directory', outfolder, ...
                                'NameValuePairs', {pchnames_k, pchvalues_k});
    sREsp3F_full{k}     = construct_fullfilename(sREsp3F, ...
                                'Directory', outfolder, ...
                                'NameValuePairs', {pchnames_k, pchvalues_k});
    sREleakF_full{k}     = construct_fullfilename(sREleakF, ...
                                'Directory', outfolder, ...
                                'NameValuePairs', {pchnames_k, pchvalues_k});
    sREparamsF_full{k}     = construct_fullfilename(sREparamsF, ...
                                'Directory', outfolder, ...
                                'NameValuePairs', {pchnames_k, pchvalues_k});
end

%% Find the IDs of cells that are stimulated or artificially activated
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

%% Build simulation commands to be read by NEURON through the here-document
sim_cmd = cell(1, ntrials);            % stores simulation commands
scmdsF_full = cell(1, ntrials);
for k = 1:ntrials
    % Update parameters for this trial
    for i = 1:length(paramvals{k})
        eval(sprintf('%s = %g;', paramnames{i}, paramvals{k}(i)));
    end

    % Commands to create neurons and build network
    sim_cmd{k} = sprintf(['buildnet("%s", %g, %g, %g, %g, %d, ', ...
                '%g, %g, %g, %g, %g, %g, %g, ', ...
                '%g, %g, %g, %g, %g, %g, %g, ', ...
                '%d, %d, %d, %d, %d, %d, %d)\n'], ...
        sREREsynF_full{k}, sp1cellID, sp2cellID, sp3cellID, ...
        REnsegs, REuseca, ...
        REcldnum, REconsyn, REtauKCC2, REepas, REdiam, REGgaba, RERErad, ...
        sp_thr, syn_del, syn_w, cai0, cao0, cli0, clo0, ...
        actcellID, actmode, savenetwork, savespikes, ...
        savesomavoltage, savesomacli, savespecial);

    % Commands to randomize leak current properties
    %     uniformly randomizes leak conductance in [REgpasLB, REgpasUB]
    sim_cmd{k} = [sim_cmd{k}, sprintf('randleak(%g, %g, "%s")\n', ...
        REgpasLB, REgpasUB, sREleakF_full{k})];

    % Commands to initialize variables
    sim_cmd{k} = [sim_cmd{k}, sprintf('vinit(%g)\n', REepas)];

    % Commands to set up neural activation protocol
    switch actmode
    case 1
        % Activate a single RE cell by injecting a train of current pulses
        sim_cmd{k} = [sim_cmd{k}, ...
                    sprintf('REsinglecp(%g, %g, %g, %g, %g, %g)\n', ...
                    actcellID, stim_start, ...
                    cp_dur, cp_amp, cp_per, cp_num)];
    case 2
        % Activate every (RERErad + 1)th RE cell by 
        %   injecting trains of current pulses
        sim_cmd{k} = [sim_cmd{k}, ...
                    sprintf('REmultcp(%g, %g, %g, %g, %g, %g, %g)\n', ...
                    actcellID, stim_start, ...
                    cp_dur, cp_amp, cp_per, cp_num, RERErad)];
    case 3
        % Activate 3 RE cells by injecting trains of current pulses
        sim_cmd{k} = [sim_cmd{k}, ...
                    sprintf('REthreecp(%g, %g, %g, %g, %g, %g, %g)\n', ...
                    actcellID, stim_start, ...
                    cp_dur, cp_amp, cp_per, cp_num, RERErad)];
    case 4
        % Activate a single RE cell at a specific voltage
        sim_cmd{k} = [sim_cmd{k}, ...
                    sprintf('REsingleact(%g, %g)\n', ...
                    actcellID, actcellv)];    
    case 5
        % Activate RE cells with a Gaussian likelihood at a specific voltage
        sim_cmd{k} = [sim_cmd{k}, ...
                    sprintf('RErandact(%g, %g, %g, %g)\n', ...
                    actcellID, actwidth, actmaxp, actcellv)];
    end

    % Commands to run simulation
    %%%%%%
    %%%%%%%%%%%%
    sim_cmd{k} = [sim_cmd{k}, sprintf(['sim(%g, %g, "%s", "%s", "%s", ', ...
                                            '"%s", "%s", "%s", ', ...
                                            '%d, %d, %d, %d)\n'], ...
        tstop, dt, ...
        sREspikeF_full{k}, sREvF_full{k}, sREcliF_full{k}, ...
        sREsp1F_full{k}, sREsp2F_full{k}, sREsp3F_full{k}, ...
        savespikes, savesomavoltage, savesomacli, savespecial)];
    %%%%%%%%%%%%
    %%%%%%

    % Commands to print all parameters
    sim_cmd{k} = [sim_cmd{k}, ...
            sprintf('print_params("%s", %g, %g, %g, %g)\n', ...
            sREparamsF_full{k}, REnsegs, REuseca, REcldnum, REconsyn)];    

    % Print simulation commands to a text file
    pchnames_k = pchnames{k};
    if iscell(pchvalues)
        pchvalues_k = pchvalues{k};
    elseif isnumeric(pchvalues)
        pchvalues_k = pchvalues(k);
    end
    scmdsF_full{k} = construct_fullfilename(scmdsF, 'Directory', outfolder, ...
                    'NameValuePairs', {pchnames_k, pchvalues_k});
    fid = fopen(scmdsF_full{k}, 'w');
    fprintf(fid, '%s\n\n', sim_cmd{k});
    fclose(fid);
end

%% Launch NEURON and execute run.hoc
results = cell(1, ntrials);         % stores simulation standard outputs
timer1 = tic();

% Get current parallel pool object without creating a new one
poolobj = gcp('nocreate');

% Create a default parallel pool object if not already exist
if isempty(poolobj)
    poolobj = parpool;
end

% Get the number of workers in the current parallel pool object
oldnumworkers = poolobj.NumWorkers;    

% Decide on the number of workers to use for running NEURON
numworkers = min(oldnumworkers, maxnumworkers_NEURON);
if renewparpoolflag_NEURON
    % Delete the parallel pool object to release memory
    delete(poolobj);
end

ct = 0;                             % counts number of trials completed
while ct < ntrials                  % while not trials are completed yet
    if singletrialnum               % if running only one trial
        % Run that trial
        first = singletrialnum;
    else
        % First trial in this batch
        first = ct + 1;             
    end
    if singletrialnum               % if running only one trial
        % Run only that trial
        last = singletrialnum;
    elseif renewparpoolflag_NEURON && ...
        ct + numworkers <= ntrials  % if memory is to be released
        % Limit the batch to numworkers
        last = ct + numworkers;
    else
        % Run all trials at once
        last = ntrials;
    end
    if renewparpoolflag_NEURON
        % Recreate a parallel pool object using 
        %   fewer workers to prevent running out of memory
        poolobj = parpool('local', numworkers);
    end
    parfor k = first:last
    %for k = first:last
        %% Use run.hoc with sim_cmd
        %##########
        %##############
        if ncells == 100
            [status, results{k}] = ...
                unix(sprintf(['x86_64/special run.hoc - << here\n', ...
                                '%s\nprint "No_Errors!"\nhere'], ...
                    sim_cmd{k}));    
        elseif ncells == 10
            [status, results{k}] = ...
                unix(sprintf(['x86_64/special run_small.hoc - << here\n', ...
                                '%s\nprint "No_Errors!"\nhere'], ...
                    sim_cmd{k}));    
        else
            status = -3;
            results{k} = 'NEURON wasn''t run because ncells is not correct\n';
        end

        fid = fopen(fullfile(outfolder, ...
                    strrep(soutF, 'sim', ['sim_', num2str(k)])), 'w');
        fprintf(fid, ['Return status was: %d\n\n', ...
                        'Simulation output was:\n\n%s\n'], status, results{k});
        fclose(fid);
        fprintf('Simulation #%d complete!\n', k);
        %##############
        %##########
   end
    if renewparpoolflag_NEURON
        % Delete the parallel pool object to release memory
        delete(poolobj);
    end
    if singletrialnum
        % Don't run again
        ct = ntrials + 1;
    else
        % Update number of trials completed
        ct = last;
    end
end
if renewparpoolflag_NEURON
    % Recreate a parallel pool object using the previous number of workers
    poolobj = parpool('local', oldnumworkers);
end

%% Analyze simulation standard outputs
time_taken = toc(timer1);
fprintf('It took %3.3g seconds to run all %d simulations with NEURON!!\n', ...
            time_taken, ntrials);
fprintf('\n');
ran_into_errors = cellfun(@isempty, strfind(results, 'No_Errors!'));
if sum(ran_into_errors) == 0
    fprintf('No Errors in NEURON!\n');
    fprintf('\n');

    % Could use the following to display simulation outputs for 
    %     debugging purposes (these are always saved as text files though)
    %{
    for k = 1:ntrials
        disp(results{k});
        fprintf('\n');
    end
    %}
else
    ind_problematic = find(ran_into_errors > 0, 1);
    fprintf(['Simulation for Sweep #', num2str(ind_problematic), ...
                ' ran into errors with output:\n']);
    error(results{ind_problematic});
    fprintf('%s\n', results{ind_problematic});
end

%% Plot stuff
timer2 = tic();
% Show network topology
if plotnetwork
    [RERE] = show_RTnet(infolder, outfolder);
end

% Show spike raster plot for each set of neurons (each .spi file)
[~, numactive, latency, actvel] ...
    = raster_plot(infolder, 'OutFolder', outfolder, ...
            'RenewParpool', renewparpoolflag_plots, ...
            'MaxNumWorkers', maxnumworkers_plots, ...
            'PlotSpikes', plotspikes, 'PlotTuning', plottuning);

% Show single neuron traces and heat maps for selected neurons 
%   (for each .singv, .singcli & .singsp file)
if plotsingleneurondata
    single_neuron(infolder, 'OutFolder', outfolder, ...
            'CellsToPlot', IDs, 'PropertiesToPlot', properties, ...
            'RenewParpool', renewparpoolflag_plots, ...
            'MaxNumWorkers', maxnumworkers_plots);
end

%% Save all variables in a mat file named by the date & time
save(fullfile(outfolder, sprintf('%s.mat', outfoldername)), '-v7.3');

%% Compute time taken
time_taken = toc(timer2);
fprintf('It took %3.3g seconds to plot and save stuff!!\n', time_taken);
fprintf('\n');

%% Play Handel if not on Rivanna
if exist('/media/adamX/RTCl/', 'dir') == 7
    load handel
    sound(y, Fs);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}