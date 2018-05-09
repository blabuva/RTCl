%% combine_curves.m
%% Combine tuning curves for the same cli0 conditions

%% Set input directories
infolders{1, 1} = '/media/adamX/RTCl/20170509T1833_cli0_5_stimfreq_1to38';
infolders{1, 2} = '/media/adamX/RTCl/20170510T0147_cli0_5_stimfreq_10to30';
infolders{2, 1} = '/media/adamX/RTCl/20170509T1836_cli0_8_stimfreq_1to38';
infolders{2, 2} = '/media/adamX/RTCl/20170510T0200_cli0_8_stimfreq_10to30';
infolders{3, 1} = '/media/adamX/RTCl/20170509T1838_cli0_11_stimfreq_1to38';
infolders{3, 2} = '/media/adamX/RTCl/20170510T0149_cli0_11_stimfreq_10to30';
infolders{4, 1} = '/media/adamX/RTCl/20170509T1840_cli0_17_stimfreq_1to38';
infolders{4, 2} = '/media/adamX/RTCl/20170510T0201_cli0_17_stimfreq_10to30';

%% Set and/or make combined folders 
comfolders{1} = '/media/adamX/RTCl/REgabaaTfall_50_REGgaba_0.0025_REtauKCC2_32_cli0_5_stimfreq_1to38';
comfolders{2} = '/media/adamX/RTCl/REgabaaTfall_50_REGgaba_0.0025_REtauKCC2_32_cli0_8_stimfreq_1to38';
comfolders{3} = '/media/adamX/RTCl/REgabaaTfall_50_REGgaba_0.0025_REtauKCC2_32_cli0_11_stimfreq_1to38';
comfolders{4} = '/media/adamX/RTCl/REgabaaTfall_50_REGgaba_0.0025_REtauKCC2_32_cli0_17_stimfreq_1to38';

%% Check if needed output directories exist
for k = 1:numel(comfolders)
    if exist(comfolders{k}, 'dir') ~= 7
        mkdir(comfolders{k});
        fprintf('New directory made: %s\n\n', comfolders{k});
    end
end

%% Copy .spi & .csv files over
for i = 1:size(infolders, 1)
    for j = 1:size(infolders, 2)
        copyfile(fullfile(infolders{i, j}, '*.spi'), comfolders{i});
        copyfile(fullfile(infolders{i, j}, '*.csv'), comfolders{i});
    end
end

%% Combine loopedparams
for i = 1:size(infolders, 1)
    combine_loopedparams(infolders{i, 1}, infolders{i, 2}, comfolders{i});
end

%% Plot tuning curves
% Output will be in /tuning_curves_across_conditions/
[pvalues, numactive, latency, actvel, nump, ...
    pnames, plabels, pislog, nperp, ncells, actmode, latency_cells_to_plot] = ...
    tuning_curves (comfolders, [], [], [], [], [], [], [], [], [], [], [], [], [], ...
            'PCondNames', {'cli0'}, ...
            'CellsToPlot', 40, ...
            'FigTypes', {'png', 'fig'});
