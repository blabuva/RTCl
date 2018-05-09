function [nump, pnames, plabels, pislog, pvalues, nperp, pchnames, pchvalues, ncells, actmode, loopmode] = get_loopedparams (infolder)
%% Get parameters that were looped in the simulation from loopedparams.mat
% USAGE: [nump, pnames, plabels, pislog, pvalues, nperp, pchnames, pchvalues, ncells, actmode, loopmode] = get_loopedparams (infolder)
% Arguments:
%     TODO
%
% Requires:
%        infolder/*loopedparams.mat
%        cd/all_ordered_pairs.m
%
% Used by:
%        cd/combine_loopparams.m
%        /media/adamX/RTCl/raster_plot.m
%        /media/adamX/RTCl/tuning_curves.m
%
% 2017-04-14 Moved from raster_plot.m
% 2017-04-17 Now extracts pvalues if it exists
% 2017-04-17 Removed latency_cells_to_plot
% 2017-05-03 Moved to Adams_Functions
% 2017-05-03 Renamed function get_loopparams -> get_loopedparams
% 2017-05-03 Changed loopedparamsfile 'loopvariables.mat' -> 'loopedparams.mat'
% 2017-05-04 Added loopedparamsfile_old for compatibility with old data

%% Must be consistent with make_loopedparams.m & combine_loopparams.m
loopedparamsfile = 'loopedparams.mat';
loopedparamsfile_old = 'loopvariables.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Find the loopedparams.mat in infolder
loopfiles = dirr(infolder, loopedparamsfile);
loopfiles_old = dirr(infolder, loopedparamsfile_old);
all_loopfiles = [loopfiles; loopfiles_old];
if length(all_loopfiles) > 1
    error('Too many versions of %s or %s!', loopedparamsfile, loopedparamsfile_old);
elseif length(all_loopfiles) < 1
    error('%s or %s is required!', loopedparamsfile, loopedparamsfile_old);
else
    % There should be only one loopedparams.mat or loopvariables.mat
    m = matfile(fullfile(infolder, all_loopfiles(1).name), ...  
        'Writable', true);          % make writable for pvalues and nperp
end

%% Extract pnames, plabels, pislog, ncells, actmode, loopmode
pnames = m.pnames;
plabels = m.plabels;
pislog = m.pislog;
ncells = m.ncells;
if isempty(whos(m, 'actmode'))      % for compatibility with older versions of loopedparams.mat
    actmode = 1;
else
    actmode = m.actmode;
end
if isempty(whos(m, 'loopmode'))     % for compatibility with older versions of loopedparams.mat
    loopmode = 'cross';
else
    loopmode = m.loopmode;
end

%% Compute nump
nump = numel(pnames);               % number of parameters changed

%% Extract pvalues if it exists and compute nperp
% Otherwise, generate parameter values used from pmin, pmax, pinc
if ~isempty(whos(m, 'pvalues'))    % pvalues already exists
    % Extract pvalues and compute nperp
    pvalues = m.pvalues;
    nperp = cellfun(@length, pvalues);

    % Check if number of parameter values is consistent
    if ~isempty(whos(m, 'nperp'))
        if ~isequal(nperp, m.nperp)
            error('Number of parameter values is inconsistent!');
        end
    elseif ~isempty(whos(m, 'ntrperp'))     % older version of loopvariables.m
        if ~isequal(nperp, m.ntrperp)
            error('Number of parameter values is inconsistent!');
        end
    end
else                % pvalues don't exist yet
    % Extract pmin, pmax & pinc
    pmin = m.pmin;
    pmax = m.pmax;
    pinc = m.pinc;

    % Generate parameter values used and compute nperp
    pvalues = cell(1, nump);
    nperp = zeros(1, nump);
    for p = 1:nump
        
        if pislog(p)
            pvalues{p} = exp(log(pmin(p)):log(pinc(p)):log(pmax(p)))';    % parameters used as a column vector
        else
            pvalues{p} = (pmin(p):pinc(p):pmax(p))';            % parameters used as a column vector
        end
        nperp(p) = length(pvalues{p});
    end

    % Save pvalues and nperp
    m.pvalues = pvalues;
    m.nperp = nperp;
end

%% Store parameter names and values for each file (to keep things in order)
if ~isempty(whos(m, 'pchnames')) && ...
    ~isempty(whos(m, 'pchvalues'))    % pchnames & pchvalues already exist
    % Extract pchnames & pchvalues
    pchnames = m.pchnames;          % changed parameter name for each file
    pchvalues = m.pchvalues;        % changed parameter value for each file
else
    switch loopmode
    case 'cross'
        % Create pchnames & pchvalues
        pchnames = cell(1, nfiles);     % changed parameter name for each file
        pchvalues = zeros(1, nfiles);   % changed parameter value for each file
        ct = 0;                         % counts number of files used
        for p = 1:nump
            % Store parameter names and values for each file (to keep things in order)
            indices = (ct + 1):(ct + nperp(p));     % indices for this parameter
            pchnames(indices) = repmat(pnames(p), 1, length(indices));
            pchvalues(indices) = pvalues{p}';

            % Update count of number of files used
            ct = ct + nperp(p);
        end
    case 'grid'
        % Repeat each parameter name nperp times
        pnames_rep = cell(1, nump);     % stores each parameter name repeated nperp times
        for p = 1:nump
            pnames_rep{p} = repmat(pnames(p), nperp(p), 1);
        end

        % Construct all possible ordered pairs of values and corresponding names
        pchvalues = all_ordered_pairs(pvalues);     % a cell array of all ordered pairs
        pchnames = all_ordered_pairs(pnames_rep);   % a cell array of ordered pairs of corresponding parameter names
    end    

    % Save pchnames & pchvalues
    m.pchnames = pchnames;
    m.pchvalues = pchvalues;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

