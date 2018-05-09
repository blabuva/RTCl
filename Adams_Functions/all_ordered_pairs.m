function [orderedPairs, nOrderedPairs, nVectors] = all_ordered_pairs (vectors, varargin)
%% Generates a cell array of all ordered pairs of elements/indices, one from each vector
% Usage: [orderedPairs, nOrderedPairs, nVectors] = all_ordered_pairs (vectors, varargin)
% Outputs:
%       orderedPairs    - a cell array of all ordered pairs of elements/indices,
%                           one from each vector(default: elements)
%       nOrderedPairs   - total number of unique ordered pairs 
%                           generated from the vectors
%       nVectors        - total number of vectors
% Arguments:    
%       vectors         - a cell array of vectors to be paired
%                       must be a cell array of vectors
%       varargin        - 'OutIndices': whether the ordered pairs are 
%                           indices instead of elements
%                       must be logical 1 (true) or 0 (false)
%                       default == false
%                       - 'VectorChangeOrder': order of vectors to change
%                       must be an unambiguous, case-insensitive match 
%                           to one of the following: 
%                           'forward'   - change first vector first, 
%                                           same as linear indicing in Matlab
%                           'backward'  - change last vector first
%                       default == 'forward'
%
% Used by:    
%       cd/make_loopedparams.m
%       cd/get_loopedparams.m
%
% File History:
% 2017-05-03 Created by Adam Lu
% 2017-05-03 Added the parameter 'VectorChangeOrder', make 'LinearIndex' default
% 2018-05-08 Changed tabs to spaces and limited width to 80
% 

%% Hard-coded parameters
possibleVectorChangeOrders = {'forward', 'backward'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to an input Parser
addRequired(iP, 'vectors', ...              % a cell array of vectors
    @(x) assert(iscell(x) && min(cellfun(@isvector, x)), ...
        'vectors must be a cell array of vectors!'));

% Add parameter-value pairs to the input Parser
addParameter(iP, 'OutIndices', false, ... 
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'VectorChangeOrder', 'forward', ...
    @(x) any(validatestring(x, possibleVectorChangeOrders)));

% Read from the input Parser
parse(iP, vectors, varargin{:});
outindices = iP.Results.OutIndices;
vecorder = validatestring(iP.Results.VectorChangeOrder, ...
                            possibleVectorChangeOrders);

%% Do the job
% Find numbers
nVectors = numel(vectors);              % number of vectors
nElements = cellfun(@length, vectors);  % number of elements in each vector
nOrderedPairs = prod(nElements);        % final number of ordered pairs

% Construct partial products
%   Note: the first/last partial product is always 1
partialProducts = ones(1, nVectors);
if nVectors > 1             % only do this if there are at least two vectors
    switch vecorder         % order of vectors to change
    case 'forward'
        % Change first vector first, same as linear indicing in Matlab
        %   Start from the second partial product
        for v = 2:nVectors
            % Each partial product is the previous partial product multiplied 
            %   by the number of elements in the previous vector
            partialProducts(v) = partialProducts(v - 1) * nElements(v - 1);
        end
    case 'backward'
        % Change last vector first
        %   Start from the next to last partial product
        for v = (nVectors - 1):-1:1
            % Each partial product is the next partial product multiplied 
            %   by the number of elements in the next vector
            partialProducts(v) = partialProducts(v + 1) * nElements(v + 1);
        end
    end
end

%% Construct ordered pairs
orderedPairsIndices = cell(1, nOrderedPairs);
                        % stores ordered pairs of indices, one from each vector
if ~outindices
    orderedPairsElements = cell(1, nOrderedPairs);
                        % stores ordered pairs of elements, one from each vector
end
for o = 1:nOrderedPairs    % for each ordered pair
    % Initialize indices and elements
    orderedPairsIndices{o} = zeros(1, nVectors);
    if isnumeric(vectors{1})
        orderedPairsElements{o} = zeros(1, nVectors);
    elseif iscell(vectors{1})
        orderedPairsElements{o} = cell(1, nVectors);
    elseif isstruct(vectors{1})
        orderedPairsElements{o} = repmat(struct, 1, nVectors);
    end

    % Find the index and element for each vector
    for v = 1:nVectors        % for each vector
        % Find the index for this vector
        orderedPairsIndices{o}(v) = ...
            mod(floor((o - 1)/partialProducts(v)), nElements(v)) + 1;

        % Find the corresponding element from this vector
        if ~outindices
            if isnumeric(vectors{1})
                orderedPairsElements{o}(v) = ...
                    vectors{v}(orderedPairsIndices{o}(v));
            elseif iscell(vectors{1}) || isstruct(vectors{1})
                orderedPairsElements{o}{v} = ...
                    vectors{v}{orderedPairsIndices{o}(v)};
            end
        end
    end
end

%% Choose orderer pair for output
if outindices
    orderedPairs = orderedPairsIndices;
else
    orderedPairs = orderedPairsElements;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}
