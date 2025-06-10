function [nn_red, S] = computeReducedNetwork(obj, S, varargin)
% computeReducedNetwork - computes a reduces network by merging similar
%   neurons based on the given input. Note that we transform the network
%   into normal form (alternating linear & nonlinear layer) prior to the
%   reduction.
%
% Syntax:
%   [nn_red, S] = neuralNetwork.computeReducedNetwork(obj, S)
%   [nn_red, S] = neuralNetwork.computeReducedNetwork(obj, S, varargin)
%
% Inputs:
%    obj - neuralNetwork
%    S - zonotope or polyZonotope
%    varargin - name-value pairs
%       <'BucketType',type> - 'static' or 'dynamic'
%       <'ReductionRate',reductionRate> - max. rate of remaining neurons
%       <'BucketTol',tol> - (initial) bucket tolerance
%       <'InputCompression',doCompression> - whether input set should be compressed 
%       <'Verbose',verbose> - verbose output
%       <'Plot',doPlot> - whether information should be plotted
%       <'MinSearches',minSearches> - min searches of binary search algorithm
%
% Outputs:
%    nn_red - reduced neural network
%    pZ - output set
%
% Reference:
%    [1] Ladner et al., Fully Automatic Neural Network Reduction for 
%        Formal Verification. arxiv. 2024.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork/getNormalForm

% Authors:       Tobias Ladner
% Written:       20-December-2022
% Last update:   09-January-2023
%                09-November-2023 (TL, CORA v2024 changes)
%                24-November-2023 (TL, added doCompression)
%                29-November-2023 (TL, major restructuring, added reductionRate)
%                29-July-2024 (TL, limit number of (dynamics) buckets)
%                27-August-2024 (TL, simplified bounds handling)
%                09-November-2024 (TL, store indices of merged neurons)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
narginchk(2,inf)
[nn,S,verbose,type,tol,doCompression,doPlot,reductionRate,minSearches] = aux_parseInput(obj,S,varargin{:});

if ~isempty(reductionRate) && min(reductionRate) == 1
    % reductionRate = 1 means no reduction; returning original network
    if verbose
        disp('No reduction. Returning original network.')
    end
    nn_red = nn;
    if nargout == 2
        S = nn.evaluate(S);
    end
    return
end

% compress input
if doCompression
    % select first entry of parameters for compression
    tol_comp = tol(1);
    reductionRate_comp = [];
    if ~isempty(reductionRate)
        reductionRate_comp = reductionRate(1);
    end

    % do input set compression
    [nn,S] = aux_compressInput(nn,S,type,tol_comp,reductionRate_comp,verbose);
    % enlarge vectors to match new network
    tol = tol([1,1,1:length(tol)]);
    if ~isempty(reductionRate)
        reductionRate = reductionRate([1,1,1:length(reductionRate)]);
    end
end

% save number of neurons of original network
numNeurons = nn.getNumNeurons();

% parameters 
options = aux_initOptions();

% compute bounds of input
bounds = interval(S);

% assuming alternating linear and nonlinear layers
kappa = length(nn.layers);
for k = 2:2:(kappa-1)
    % extract layers
    nnLinIn = nn.layers{k-1};
    nnAct = nn.layers{k};
    nnLinOut = nn.layers{k+1};
    
    % propagate bounds to next activation layer k
    options.nn.reuse_bounds = false;
    bounds_pre = nn.evaluate(bounds, options, k-1:k);
    bounds_pre_inf = bounds_pre.inf;
    bounds_pre_sup = bounds_pre.sup;
    
    if doPlot
        figure; hold on;
        histogram(bounds_pre.inf, 100)
        histogram(bounds_pre.sup, 100)
    end

    % determine buckets
    buckets = [];
    if strcmp(type, 'static')
        buckets = nnAct.getMergeBuckets();
    elseif strcmp(type, 'dynamic')
        buckets = center(bounds_pre)';
        index = randperm(numel(buckets));
        buckets = buckets(index);
    end
        
    % limit number of buckets
    buckets = buckets(1:min(100,numel(buckets)));

    % init used parameters for layer k
    tol_k = tol(k);
    tol_k_min = 0; tol_k_max = max(bounds_pre_sup)-min(bounds_pre_inf);
    if isempty(tol_k_max)
        throw(CORAerror('CORA:specialError','Invalid upper bounds for binary search.'))
    end
    minSearches_k = minSearches;
    while true % reduce network until desired rate is reached
        minSearches_k = minSearches_k - 1;

        % bucket bounds
        bInf = buckets-tol_k;
        bSup = buckets+tol_k;

        % compute containment and filter only first belonging
        M_merged = (bInf <= bounds_pre_inf) & (bounds_pre_sup <= bSup);
        M_merged = (cumsum(M_merged, 2) == 1) & (M_merged == 1);
        
        % only select buckets with more than one containments
        idx_b = sum(M_merged, 1) > 1;
        M_merged = M_merged(:, idx_b);
        
        % store results
        M_merged = M_merged';

        % mark chosen neurons
        idx = any(M_merged, 1);
        
        if isempty(reductionRate)
            % no desired rate specified
            break;
        else
            % compute desired and actual reduction rate
            rate_k_desired = reductionRate(k);
            rate_k = sum(~idx) / length(idx);
            
            % check if actual rate is ok
            if rate_k_desired == rate_k || ...
                (minSearches_k <= 0 && rate_k <= rate_k_desired)
                % reduction done
                break;
            end

            % scale bucket tolerance to get closer to desired rate
            % via binary search
            if rate_k < rate_k_desired
                % too much merged -> reduce bucket tolerance
                tol_k_max = tol_k;
                tol_k = (tol_k_min + tol_k)/2;
            else
                % not enough merged -> increase bucket tolerance
                tol_k_min = tol_k;
                if ~isinf(tol_k_max)
                    tol_k = (tol_k + tol_k_max)/2;
                else
                    % rapidly enlarge bucket tolerance 
                    % to quickly reach too large bucket tolerance
                    tol_k = 10*tol_k;
                end
            end

            % stop after min and max tolerance converged
            % (usually only happens for compressed input)
            if withinTol(tol_k_min,tol_k_max,1e-10)
                % choose last tolerance that was too large
                tol_k = tol_k_max;
            end
        end
    end
    % update tolerance vector with tolerance used
    tol(k) = tol_k;

    % count number of merged neurons
    num_merged = size(M_merged, 1);

    % init merge matrix
    M_unmerged = diag(sparse(~idx)); % keep un-merged neurons
    M_unmerged = M_unmerged(any(M_unmerged, 2), :); % delete zero rows
    
    % merge 'input' weight matrix
    W1 = nnLinIn.W;
    b1 = nnLinIn.b;
    W1m = M_unmerged * W1;
    b1m = M_unmerged * b1;
    
    % merge 'output' weight matrix
    W2 = nnLinOut.W;
    b2 = nnLinOut.b;
    W2m = W2 * M_unmerged'; % sum
    b2m = b2;  % unchanged!

    % init linear layers \widehat{L}_{k-1}, \widehat{L}_{k+1}
    nnLinInNew = nnLinearLayer(full(W1m), full(b1m), sprintf('%s_reduced',nnLinIn.name));
    nnActNew = nnAct.copy();
    nnActNew.name = [nnActNew.name(1:end-5) '_reduced'];
    nnLinOutNew = nnLinearLayer(full(W2m), b2m, sprintf('%s_reduced',nnLinOut.name));

    % compute approx error
    if num_merged > 0
        % select bounds of merged neurons
        approx_error = bounds_pre .* idx';
        % propagate forward
        nnLinOutNew.d = W2*approx_error;

        % uncomment this line if you want to use 
        % the value from the merge bucket:
        % nnLinOutNew.d = sum(sum(((cumsum(M_merged')' == M_merged) & (M_merged == 1)) .* approx_error',2) .* M_merged,1)';
    end

    if ~representsa_(nnLinOut.d,'emptySet',eps)
        % add approx error from previous reduction
        if representsa_(nnLinOutNew.d,'emptySet',eps)
            nnLinOutNew.d = nnLinOut.d;
        else
            nnLinOutNew.d = nnLinOutNew.d + nnLinOut.d;
        end
    end

    % keep old approx error in L_{k+1} for unmerged dimensions
    if ~representsa_(nnLinIn.d,'emptySet',eps)
        nnLinInNew.d = nnLinIn.d(~idx);
    end

    % store which neurons were merged in new activation layer
    nnActNew.merged_neurons = idx;
     
    % update layers
    nn.layers{k-1} = nnLinInNew;     % \widehat{L}_{k-1}
    nn.layers{k} = nnActNew;         % \widehat{L}_{k}
    nn.layers{k+1} = nnLinOutNew;    % \widehat{L}_{k+1}
    
    % evaluate set on reduced network
    S = nn.evaluate(S, options, k-1:k);

    % update bounds
    bounds = interval(S);
end

% propagate through output layers
S = nn.evaluate(S, options, k+1:kappa);

% compute reduction rate
numNeuronsRed = nn.getNumNeurons();
rate = sum(numNeuronsRed(2:2:end-1))/sum(numNeurons(2:2:end-1));

if isempty(reductionRate)
    nn.reductionRate = numNeuronsRed ./ numNeurons;
else
    nn.reductionRate = reductionRate;    
end

if verbose
    % display resulting number of neurons
    disp([numNeurons; numNeuronsRed])
    fprintf("Remaining neurons within network: %.2f%%\n", rate*100)
end

% sanity check
% N = 500;
% xs = pZ.randPoint(N);
% ys = nn.evaluate(xs);
% bounds = nn.evaluate(bounds, options, 13);
% 
% res = false(1, N);
% for i=1:N
%     res(i) = bounds.contains(ys(:, i));
% end
% disp(all(res))

if doPlot
    figure; hold on;
    histogram(bounds.inf, 100)
    histogram(bounds.sup, 100)
end

nn_red = nn;

end


% Auxiliary functions -----------------------------------------------------

function [nn,S,verbose,type,tol,doCompression,doPlot,reductionRate,minSearches] = aux_parseInput(obj,S,varargin)

% name-value pairs -> number of input arguments is always a multiple of 2
if mod(nargin,2) ~= 0
    throw(CORAerror('CORA:evenNumberInputArgs'));
else
    % read input arguments
    NVpairs = varargin;
    % check list of name-value pairs
    checkNameValuePairs(NVpairs,{'Verbose','BucketType','BucketTol','InputCompression','Plot','ReductionRate','MinSearches'});
    [NVpairs,verbose] = readNameValuePair(NVpairs,'Verbose',{'islogical','isscalar'},false);
    [NVpairs,type] = readNameValuePair(NVpairs,'BucketType',{},'static');
    [NVpairs,tol] = readNameValuePair(NVpairs,'BucketTol',{'isnumeric','isscalar'},1e-2);
    [NVpairs,doCompression] = readNameValuePair(NVpairs,'InputCompression',{'islogical','isscalar'},false);
    [NVpairs,doPlot] = readNameValuePair(NVpairs,'Plot',{'islogical','isscalar'},false);
    [NVpairs,reductionRate] = readNameValuePair(NVpairs,'ReductionRate',{'isnumeric'},[]);
    [NVpairs,minSearches] = readNameValuePair(NVpairs,'MinSearches',{'isnumeric','isscalar'},10);
end

% bring network in normal form (and create copy)
nn = obj.getNormalForm();
kappa = length(nn.layers);

% check bucket type
methods = {'static','dynamic'};
if ~ismember(type, methods)
    throw(CORAerror("CORA:wrongValue","name-value pair 'BucketType'",methods))
end

% reduction rate
if isscalar(reductionRate)
    % extend to number of layers
    reductionRate = reductionRate * ones(1,kappa);
end
if ~isempty(reductionRate) && length(reductionRate) ~= kappa
   throw(CORAerror('CORA:wrongValue',"name-value pair 'ReductionRate'",'scalar or length equal to number of layers')) 
end

% bucket tolerance
if isscalar(tol)
    % extend to number of layers
    tol = tol * ones(1,kappa);
end
if ~isempty(tol) && length(tol) ~= kappa
   throw(CORAerror('CORA:wrongValue',"name-value pair 'BucketTol'",'empty, scalar, or length equal to number of layers')) 
end


end

function [nn,S] = aux_compressInput(nn,S,type,tol,rate,verbose)

% input has to be positive for the following code to work
lb = interval(S).inf;
if any(lb < 0)
    throw(CORAerror("CORA:specialError",'Input has to be positive for input compression.'))
end


 % construct identity network
neurons_in = dim(S);
nnComp = neuralNetwork({ ...
    nnLinearLayer(eye(neurons_in)); ...
    nnReLULayer();
    nnLinearLayer(eye(neurons_in)); ...
});

if verbose
    disp('Compressing input ...')
end

% reduce
nnComp = nnComp.computeReducedNetwork(S,"BucketType",type,"BucketTol",tol,"ReductionRate",rate);

% construct new input
S = nnComp.evaluate(S,struct,1);
S = zonotope(S);

nn = neuralNetwork([ ...
    nnComp.layers(3); ...
    {nnReLULayer()}; ... % necessary to keep alternating layers
    nn.layers
]);

end

function options = aux_initOptions()

options = struct();
options.nn.poly_method = "regression";
options.nn.bound_approx = true;
% options.nn.num_generators = 50000;
options.nn.propagate_bounds = false;
options.nn.do_pre_order_reduction = false;
options.nn.remove_GI = false;
options.nn.add_approx_error_to_GI = true;

end

% ------------------------------ END OF CODE ------------------------------
