function [RERE] = show_RTnet (infolder, outfolder)
%% Shows network topology for each RT network (each .syn file in the infolder)
% USAGE: [RERE] = show_RTnet (infolder, outfolder)
%
% Arguments:
%       infolder   - the name of the directory containing the .syn files, e.g. '20170317T1127_Ggaba_0.01'
%                   must be a character array
%       outfolder  - (opt) the name of the directory that the plots will be placed
%                   must be a character array
%                   default: same as infolder
% Requires:
%       infolder/*.syn
%       /home/Matlab/Downloaded_Functions/dirr.m
% Used by:
%       cd/neuronlaunch.m
%
% 2017-02-08 Modified from show_net.hoc of MCM 201505 for Peter's project
% 2017-02-13 Now creates a network plot for multiple RT networks; added dirr()
% 2017-02-13 Added dirr()
% 2017-02-23 Added a zoomed-in version of the network plot
% 2017-03-14 Changed the way ylimits are defined
% 2017-03-27 Filenames are now modified in neuronlaunch.m


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check arguments
if nargin < 1
    error('An infolder is required, type ''help show_RTnet'' for usage');
elseif isempty(infolder) || ~isdir(infolder)
    error('infolder must be a directory!');
elseif nargin >= 2 && ~isdir(outfolder)
    error('outfolder must be a directory!');
end

%% Set default arguments
if nargin < 2
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
addpath(fullfile(functionsdirectory, '/Downloaded_Functions/'));
                                                            % for dirr.m

%% Find all .syn files
files = dirr(infolder, '.syn');
nfiles = length(files);    

%% Create raster plot for each file
RERE = cell(nfiles, 1);        % raw network topology
parfor i = 1:length(files)
    % Load data
    RERE{i} = load(fullfile(infolder, files(i).name));

    % Set figure name (without .png)
    figname = fullfile(outfolder, strrep(files(i).name, '.syn', '_network_topology.png'));

    % Plot data
    network_topology(RERE{i}, figname, 'full', files(i).name);
    network_topology(RERE{i}, figname, 'part', files(i).name);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function network_topology(data, figname, plotmode, filename)
%% Plot network topology from data

% Prepare the figure
h = figure(floor(rand()*10^4));
clf(h);
hold on;

% Define color map
cm = colormap(lines);
ncolors = size(cm, 1);

% Plot each line
for i=1:length(data)
    if strcmp(plotmode, 'full')
        plot([0.25; 0.75], round([data(i, 1); data(i, 2)]), ':', 'Color', cm(mod(data(i, 1),ncolors) + 1, :));
    elseif strcmp(plotmode, 'part')
        plot([0.25; 0.75], round([data(i, 1); data(i, 2)]), '-', 'Color', cm(mod(data(i, 1),ncolors) + 1, :));
    end
end

% Restrict axes range
ymax = max(data(:, 1)) + 1;
ymin = min(data(:, 1)) - 1;
if strcmp(plotmode, 'full')
    axis([0 1 ymin ymax]);
elseif strcmp(plotmode, 'part')
    axis([0.2 0.8 ymin ymin+11]);
end

% Set labels
ax = gca;
set(ax, 'XTick', [0.25, 0.75]);
set(ax, 'XTickLabel', {'From', 'To'});
ylabel('Neuron number')
title('Network topology of RT network')
title(['Network topology of RT network for ', strrep(filename, '_', '\_')]);

% Save figure
if strcmp(plotmode, 'part')
    figname = strrep(figname, 'topology', 'topology_zoomed_in');
end
saveas(h, figname, 'png');

% Close figure
% close(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}
