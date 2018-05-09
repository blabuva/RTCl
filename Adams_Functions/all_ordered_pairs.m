function [ordered_pairs, num_ordered_pairs, num_vectors] = all_ordered_pairs (vectors, varargin)
%% Generates a cell array of all ordered pairs of elements/indices, one from each vector
% Usage: [ordered_pairs, num_ordered_pairs, num_vectors] = all_ordered_pairs (vectors, varargin)
% Outputs:	ordered_pairs		- a cell array of all ordered pairs of elements/indices, one from each vector
%						(default: elements)
%		num_ordered_pairs	- total number of unique ordered pairs generated from the vectors
%		num_vectors		- total number of vectors
% Arguments:	vectors		- a cell array of vectors to be paired
%				must be a cell array of vectors
% 		varargin	- 'OutIndices': whether the ordered pairs are indices instead of elements
%				must be logical 1 (true) or 0 (false)
%				default == false
%				- 'VectorChangeOrder': order of vectors to change
%				must be an unambiguous, case-insensitive match to one of the following: 
%					'forward'	- change first vector first, same as linear indicing in Matlab
%					'backward'	- change last vector first
%				default == 'forward'
%
% Used by:	
%		cd/make_loopedparams.m
%       cd/get_loopedparams.m
%
% File History:
% 2017-05-03 Created by Adam Lu
% 2017-05-03 Added the parameter 'VectorChangeOrder', make 'LinearIndex' default
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
	error('Not enough input arguments, type ''help all_ordered_pairs'' for usage');
end

% Add required inputs to an input Parser
iP = inputParser;
addRequired(iP, 'vectors', ...				% a cell array of vectors
	@(x) assert(iscell(x) && min(cellfun(@isvector, x)), ...
		'vectors must be a cell array of vectors!'));

% Add parameter-value pairs to the input Parser
addParameter(iP, 'OutIndices', false, ...		% whether the ordered pairs are indices instead of elements
	@(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'VectorChangeOrder', 'forward', ...	% order of vectors to change
	@(x) any(validatestring(x, {'forward', 'backward'})));

% Read from the input Parser
parse(iP, vectors, varargin{:});
outindices = iP.Results.OutIndices;
vecorder = validatestring(iP.Results.VectorChangeOrder, {'forward', 'backward'});

%% Find numbers
num_vectors = numel(vectors);			% number of vectors
num_elements = cellfun(@length, vectors);	% number of elements in each vector
num_ordered_pairs = prod(num_elements);		% final number of ordered pairs

%% Construct partial products
partial_products = ones(1, num_vectors);	% the first/last partial product is always 1
if num_vectors > 1			% only do this if there are at least two vectors
	switch vecorder			% order of vectors to change
	case 'forward'			% change first vector first, same as linear indicing in Matlab
		for v = 2:num_vectors		% start from the second partial product
			% Each partial product is the previous partial product multiplied by the number of elements in the previous vector
			partial_products(v) = partial_products(v - 1) * num_elements(v - 1);
		end
	case 'backward'			% change last vector first
		for v = (num_vectors - 1):-1:1	% start from the next to last partial product
			% Each partial product is the next partial product multiplied by the number of elements in the next vector
			partial_products(v) = partial_products(v + 1) * num_elements(v + 1);
		end
	end
end

%% Construct ordered pairs
ordered_pairs_indices = cell(1, num_ordered_pairs);		% stores ordered pairs of indices, one from each vector
if ~outindices
	ordered_pairs_elements = cell(1, num_ordered_pairs);	% stores ordered pairs of elements, one from each vector
end
for o = 1:num_ordered_pairs	% for each ordered pair
	% Initialize indices and elements
	ordered_pairs_indices{o} = zeros(1, num_vectors);
	if isnumeric(vectors{1})
		ordered_pairs_elements{o} = zeros(1, num_vectors);
	elseif iscell(vectors{1})
		ordered_pairs_elements{o} = cell(1, num_vectors);
	elseif isstruct(vectors{1})
		ordered_pairs_elements{o} = repmat(struct, 1, num_vectors);
	end

	% Find the index and element for each vector
	for v = 1:num_vectors		% for each vector
		% Find the index for this vector
		ordered_pairs_indices{o}(v) = mod(floor((o - 1)/partial_products(v)), num_elements(v)) + 1;

		% Find the corresponding element from this vector
		if ~outindices
			if isnumeric(vectors{1})
				ordered_pairs_elements{o}(v) = vectors{v}(ordered_pairs_indices{o}(v));
			elseif iscell(vectors{1}) || isstruct(vectors{1})
				ordered_pairs_elements{o}{v} = vectors{v}{ordered_pairs_indices{o}(v)};
			end
		end
	end
end

%% Choose orderer pair for output
if outindices
	ordered_pairs = ordered_pairs_indices;
else
	ordered_pairs = ordered_pairs_elements;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:
TODO: Place older versions of the code that you want to save here, in case you need it back in the future

mod(floor(k/1), nperp(3))
mod(floor(k/nperp(3)), nperp(2))
mod(floor(k/(nperp(2)*nperp(3))), nperp(1))
3 * nperp(2) * nperp(3) + 6 * nperp(3) + 1 * 1


%}
