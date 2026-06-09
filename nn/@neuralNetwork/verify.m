function [res, x_, y_] = verify(nn, x, r, A, b, safeSet, varargin)
% verify - automated verification for specification on neural networks.
%
% Syntax:
%    [res, z] = nn.verify(x, r, A, b, options)
%
% Inputs:
%    nn - object of class neuralNetwork
%    x, r - center and radius of the initial set (can already be a batch)
%    A, b - specification, prove A*y <= b
%    safeSet - bool, safe-set or unsafe-set [optional: if > 1 then the
%       number of union (safe set) constraints.]
%    options - struct, evaluation options
%    timeout - positive integer, timeout in seconds
%    verbose - bool, print verbose output
%    plotDims - 2x2 plot dimensions; empty for no plotting;
%       plotDims(1,:) for input and plotDims(2,:) for output; sets
%       are stored in res.Xs, res.uXs
%
% Outputs:
%    res - result: true if specification is satisfied, false if not, empty if unknown
%    x_ - counterexample in terms of an initial point violating the specs
%    y_ - output for x_
%
% References:
%    [1] VNN-COMP'24
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Lukas Koller
% Written:       23-November-2021
% Last update:   14-June-2024 (LK, rewritten with efficient splitting)
%                20-January-2025 (LK, constraint zonotope splitting)
%                07-April-2026 (LK, memory usage improvements)
%                02-June-2026 (BK, optional input-side polytope constraints)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Check number of input arguments.
narginchk(6,13);

% Validate parameters.
[options, timeout, verbose, plotDims, ~, A_in, b_in] = ...
    setDefaultValues({struct, 100, false, [], false, [], []}, varargin);
plotting = ~isempty(plotDims);

% Validate parameters.
inputArgsCheck({ ...
    {nn, 'att','neuralNetwork'}; ...
    {x, 'att',{'numeric','gpuArray'}}; ...
    {r, 'att',{'numeric','gpuArray'}}; ...
    {A, 'att',{'numeric','gpuArray'}}; ...
    {b, 'att',{'numeric','gpuArray'}}; ...
    {options,'att','struct'}; ...
    {timeout,'att','numeric','scalar'}; ...
    {verbose,'att','logical'}; ...
    {plotting,'att','logical'}; ... % TODO: re-implement plotting
    })
options = nnHelper.validateNNoptions(options,true);

nSplits = options.nn.num_splits; % Number of input splits per dimension.
nDims = options.nn.num_dimensions; % Number of input dimension to splits.
nNeur = options.nn.num_neuron_splits; % Number of neurons to split.
nReLU = options.nn.num_relu_constraints; % Number of ReLU constraints.

if isinf(nReLU)
    % Compute the maximum number of ReLU constraints.
    nReLU = sum(cellfun(@(li) ...
        isa(li,'nnReLULayer')*prod(li.getOutputSize(li.inputSize)), ...
        nn.layers) ...
        );
end

% Extract parameters.
bSz = options.nn.train.mini_batch_size;

% Obtain number of input dimensions.
n0 = size(x,1);
% Limit the number of dimensions to split.
nDims = min(nDims,n0);
% Check the maximum number of input generators.
numInitGens = min(options.nn.train.num_init_gens,n0);
% Obtain the maximum number of approximation errors in an activation layer.
nk_max = max(cellfun(@(li) ...
    isa(li,'nnActivationLayer')*prod(li.getOutputSize(li.inputSize)), ...
    nn.layers));
% We always have to use the approximation during set propagation to ensure
% soundness.
options.nn.use_approx_error = true;
% Ensure the interval-center flag is set, if there are less generators than
% input dimensions.
options.nn.interval_center = ...
    (options.nn.train.num_approx_err < nk_max) || (numInitGens < n0);

% Specify indices of layers for propagation.
idxLayer = 1:length(nn.layers);
% Enumerate the layers of the neural networks.
[layers,~,~,~] = nn.enumerateLayers();

% To speed up computations and reduce GPU memory, we only use single
% precision.
inputDataClass = single(1);
% Check if a GPU is used during verification.
if options.nn.train.use_gpu
    % Training data is also moved to GPU.
    inputDataClass = gpuArray(inputDataClass);
end
% Move weights of the neural network to GPU.
nn.castWeights(inputDataClass);
% Clear all stored variables in all layers.
layers = aux_clearLayerFields(layers);

% Move the specification ot the GPU.
A = cast(A,'like',inputDataClass);
b = cast(b,'like',inputDataClass);

% In each layer, store ids of active generators and identity matrices
% for fast adding of approximation errors.
q = nn.prepareForZonoBatchEval(x,options,idxLayer);

% Initialize queue.
xs = x;
rs = r;
nrXs = zeros([0 size(x,2)]);
% Compute number of union constraints (all intersection constrains ->
% union of one single constraint).
if safeSet > 1
    numUnionConst = safeSet;
elseif safeSet == 1
    numUnionConst = size(A,1);
else
    % There are no union constraints.
    numUnionConst = 1;
end
% Reset to be a bool.
safeSet = logical(safeSet);
% Initialize result.
res.str = 'UNKNOWN';
x_ = [];
y_ = [];

% Initialize iteration stats.
numVerified = 0;
% Initialize iteration counter.
iter = 1;

% Specify the batch variables so we can clear them after every iteration.
batchVars = {'xi','ri','nrXi','S_','S','sens','cxi','Gxi',...
    'yic','yid','Gyi','ld_yi','ld_Gyi'};

if verbose
    % Setup table.
    table = CORAtable('double', ...
        {'Iteration','#Queue','#Verified Branches','avg. Radius', ...
        '~Unknown Vol. [%]'}, ...
        {'d','d','d','.3e','.4f'});
    table.printHeader();
end

% Specify the heuristics (see computeHeuristics for options).
% Heuristic for selecting input generators.
inputGenHeuristic = options.nn.input_generator_heuristic;
% Heuristic for splitting input dimensions.
inputSplitHeuristic = options.nn.input_split_heuristic;
% Heuristic for neuron splitting.
neuronSplitHeuristic = options.nn.neuron_split_heuristic;
if nNeur == 0
    % We are not splitting any neurons; we do not need a heuristic.
    neuronSplitHeuristic = 'none';
end
% Heuristic for ReLU tightening constraints.
reluConstrHeuristic = options.nn.relu_constraint_heuristic;
if nReLU == 0
    % We are not splitting any neurons; we do not need a heuristic.
    reluConstrHeuristic = 'none';
end

heuristics = {inputGenHeuristic,inputSplitHeuristic,...
    neuronSplitHeuristic,reluConstrHeuristic};

% The sensitivity is used for selecting input generators, neuron
% -splitting, and FGSM attacks.
computeSens = ...
    any(strcmp('gap-sensitivity',heuristics)) || ...
    any(strcmp('product-sensitivity',heuristics)) || ...
    any(strcmp('centered-sensitivity',heuristics)) || ...
    any(strcmp('most-sensitive-approx-error',heuristics)) || ...
    any(strcmp('most-sensitive-input-radius',heuristics)) || ...
    any(strcmp('most-unstable',heuristics)) || ...
    any(strcmp('least-unstable',heuristics)) || ...
    strcmp(options.nn.falsification_method,'fgsm') || ...
    (strcmp(options.nn.approx_error_order,'sensitivity*length') ...
    && (options.nn.train.num_approx_err < nk_max));
% Store the approximation error for heuristics computations.
options.nn.store_approx_error = ...
    any(strcmp('most-sensitive-approx-error',heuristics));
% Store the approximation error gradients for heuristics computations.
computeGrads = ...
    any(strcmp('zono-norm-gradient',heuristics)) || ...
    any(strcmp('least-unstable-gradient',heuristics));
options.nn.store_approx_error_grad = ...
    any(strcmp('zono-norm-gradient',{neuronSplitHeuristic,reluConstrHeuristic})) || ...
    any(strcmp('least-unstable-gradient',{neuronSplitHeuristic,reluConstrHeuristic}));
options.nn.backprop_without_weight_update = ...
    any(strcmp('zono-norm-gradient',heuristics)) || ...
    any(strcmp('least-unstable-gradient',heuristics));

timerVal = tic;

% Main splitting loop, i.e., while the queue is non empty.
while size(xs,2) > 0

    try % Catch any error and adapt parameters accordingly.

        % Check if we reach the maximum number of iterations.
        if iter > options.nn.max_verif_iter
            break;
        end

        % Obtain the current time.
        time = toc(timerVal);
        % Check the timeout.
        if time > timeout
            % Time is up.
            res.time = time;
            break;
        end

        if verbose
            % Print the iteration stats.
            aux_printIterationStats(table,iter,xs,r,rs,numVerified);
        end

        % Pop next batch from the queue.
        [xi,ri,nrXi,xs,rs,nrXs] = aux_pop(xs,rs,nrXs,bSz,options);

        % Move the batch to the GPU.
        xi = cast(xi,'like',inputDataClass);
        ri = cast(ri,'like',inputDataClass);
        nrXi = cast(nrXi,'like',inputDataClass);

        % 1. Verification -------------------------------------------------
        % 1.1. Compute all per-batch artifacts via the helper.
        [~,~,~,~,~,ld_yi,ld_Gyi,ld_Gyi_err,~,~,~,~] = ...
            aux_computeBatchArtifacts(nn,options,layers,idxLayer,xi,ri, ...
            nrXi,A,b,q,numInitGens,inputGenHeuristic,computeSens);

        % Compute the radius of the logit difference.
        ld_ri = sum(abs(ld_Gyi),2) + ld_Gyi_err;
        % 2.3. Check specification.
        if safeSet
            % safe iff all(A*y <= b) <--> unsafe iff any(A*y > b)
            % Thus, unknown if any(A*y > b).
            unknown = any(ld_yi + ld_ri(:,:) > b,1);
        else
            % unsafe iff all(A*y <= b) <--> safe iff any(A*y > b)
            % Thus, unknown if all(A*y <= b).
            unknown = all(ld_yi - ld_ri(:,:) <= b,1);
        end

        % Update counter for verified patches.
        numVerified = numVerified + sum(~unknown,'all');

        if all(~unknown)
            % Verified all subsets of the current batch. We can skip to
            % next iteration.
            iter = iter + 1;
            continue;
        elseif any(~unknown)
            % Only keep un-verified batch entries.
            xi(:,~unknown) = [];
            ri(:,~unknown) = [];
            nrXi(:,~unknown) = [];
        end
        clearvars('ld_yi')

        % Recompute all batch artifacts for the unverified entries.
        [cxi,Gxi,inputDimIds,yi,Gyi,~,ld_Gyi,~,S,sens,grad,a] = aux_computeBatchArtifacts( ...
            nn,options,layers,idxLayer,xi,ri,nrXi,A,b,q,numInitGens,inputGenHeuristic, ...
            computeSens,computeGrads,nNeur,nReLU,nSplits);
        % Obtain the current batch size.
        [~,~,cbSz] = size(Gxi);

        % 2. Falsification ----------------------------------------------------

        % 2.1. Compute adversarial examples.
        switch options.nn.falsification_method
            case 'fgsm'
                % Obtain number of constraints.
                [p,~] = size(A);
                % Try to falsification with a FGSM attack.
                if safeSet
                    dy = pagemtimes(-A,S);
                    % We combine all constraints for a stronger attack.
                    p = 1;
                else
                    dy = pagemtimes(A,S);
                end
                % If there are multiple output constraints we try to falsify
                % each one individually.
                sdy = reshape(permute(sign(dy),[2 3 1]),[n0 cbSz*p]);

                % Compute adversarial attacks based on the sensitivity.
                xi_ = repelem(xi,1,p) + repelem(ri,1,p).*sdy;

                % Clear unused variables.
                clear('dy','sdy');
            case 'center'
                % Use the center for falsification.
                xi_ = xi;
            case 'zonotack'
                % Obtain number of constraints.
                [p,~] = size(A);
                % Compute the vertex that minimizes the distance to each
                % halfspace.
                beta_ = -permute(sign(ld_Gyi(:,1:numInitGens,:)),[2 4 1 3]);
                if safeSet
                    % We have to reverse the signs for safe sets.
                    beta_(:,:,1:numUnionConst,:) = -beta_(:,:,1:numUnionConst,:);
                end
                % Put multiple candidates into the batch.
                beta = reshape(beta_,[numInitGens 1 p*cbSz]);

                % Compute attack.
                delta = pagemtimes(repelem(Gxi(:,1:numInitGens,:),1,1,p),beta);
                % Compute candidates for falsification.
                xi_ = repelem(xi,1,p) + delta(:,:);

                % Clear unused variables.
                clear('beta_','beta','delta');
            otherwise
                % Invalid option.
                throw(CORAerror('CORA:wrongFieldValue', ...
                    'options.nn.falsification_method', ...
                    {'fgsm','center','zonotack'}));
        end

        % 2.2. Check the specification for adversarial examples.

        % Check the adversarial examples.
        [~,~,falsified,x_,y_] = ...
            aux_checkPoints(nn,options,idxLayer,A,b,safeSet,xi_);

        if any(falsified)
            % Validate CE against input-side constraints (polytope input sets).
            if ~isempty(A_in) && ~isempty(x_)
                if any(A_in * x_ > b_in + 1e-5)
                    falsified = false(size(falsified));
                    x_ = []; y_ = [];
                end
            end
        end
        if any(falsified)
            % Found a counterexample.
            res.str = 'COUNTEREXAMPLE';
            break;
        end

        % 3. Input Refinement. --------------------------------------------

        switch options.nn.refinement_method
            case 'naive'
                % The sets are not refined; split the input dimensions.
                xis = xi;
                ris = ri;
                % Store the indices of the split dimensions.
                dimIds = NaN([0 cbSz],'like',xis);
                for i=1:nDims
                    % Compute the heuristic.
                    his = nnHelper.computeHeuristic(inputSplitHeuristic, ...
                        xis - ris, ... lower bound
                        xis + ris, ... upper bound
                        ris, ... approximation error
                        sens, ... sensitivity
                        grad, ... zonotope norm gradient
                        [],[],false,1);
                    % Split the input sets along one dimensions.
                    [xis,ris,dimId] = aux_split(xis,ris,his,nSplits);
                    % Append the split dimension.
                    dimIds = [repmat(dimIds,1,nSplits); dimId];
                    % Replicate sensitivity and criticallity value.
                    sens = repmat(sens,1,nSplits);
                    grad = repmat(grad,1,nSplits);
                end
                % There is no neuron splitting.
                nrXis = zeros([0 size(xis,2)],'like',inputDataClass);
            case 'zonotack'

                % Construct neuron-split constraints.
                if nNeur > 0 && nSplits > 1
                    % Create split constraints for neurons within the
                    % network. To save memory, the constraint were already
                    % computed during the forward propagation.

                    % Extract the computed constraints.
                    An = a.As;
                    bn = a.bs;
                    newNrXi = a.nrSplitIdx;
                else
                    % There are no general-split constraints.
                    An = zeros([0 q cbSz],'like',inputDataClass);
                    bn = zeros([0 1 cbSz],'like',inputDataClass);
                    newNrXi = -ones([0 cbSz],'like',inputDataClass);
                end

                % Identify dummy splits; we insert dummy splits to maintain
                % batch size in the case where no meaningful split can be
                % found. TODO: handle these properly; remove.
                isDummySplit = all(isinf(newNrXi),1) & any(isinf(newNrXi),1);

                % Construct input split constraints.
                if nDims > 0 && nSplits > 1
                    % When not all input dimensions get an assigned generator
                    % we have to restrict and reorder the dimensions.
                    % Therefore, we compute indices.
                    permIdx = reshape(sub2ind(size(xi), ...
                        inputDimIds,repelem(1:cbSz,numInitGens,1)), ...
                        [numInitGens cbSz]);
                    % Permute the input and radius.
                    xi_ = xi(permIdx);
                    ri_ = ri(permIdx);
                    % Permute the sensitivity.
                    if ~isempty(sens)
                        sens_ = sens(permIdx);
                    else
                        % There is no sensitivity.
                        sens_ = [];
                    end

                    % Compute the heuristic.
                    hi = nnHelper.computeHeuristic(inputSplitHeuristic, ...
                        xi_ - ri_, ... lower bound
                        xi_ + ri_, ... upper bound
                        ri_, ... approximation error
                        sens_, ... sensitivity
                        grad, ... zonotope norm gradient
                        [],[],false,1);

                    % Compute input-split constraints.
                    [Ai,bi,~,~] = ...
                        aux_dimSplitConstraints(hi(:,:),nSplits,nDims);
                else
                    % There are no input-split constraints.
                    Ai = zeros([0 q cbSz],'like',inputDataClass);
                    bi = zeros([0 1 cbSz],'like',inputDataClass);
                end

                % Pad offsets if there are different number of offsets in
                % neuron split and input split constraints.
                if size(bn,2) ~= size(bi,2)
                    bn = [bn NaN([nNeur max(size(bi,2)-size(bn,2),0) cbSz],'like',bn)];
                    bi = [bi NaN([nDims max(size(bn,2)-size(bi,2),0) cbSz],'like',bi)];
                end
                % Append zeros for generators.
                An_ = [An zeros([size(An,1) q-size(An,2) cbSz],'like',An)];
                Ai_ = [Ai zeros([size(Ai,1) q-size(Ai,2) cbSz],'like',Ai)];
                % Concatenate input and neuron splits.
                As = [An_; Ai_];
                bs = [bn; bi];
                % Pad the neuron split indices with NaN for the input dimensions.
                newNrXi = [newNrXi; NaN(size(Ai,[1 3]),'like',newNrXi)];

                if nReLU > 0
                    % Construct and reduce the aggregated relu constraints.
                    [At,bt] = aux_extractAndReduceReluConstraints(a);
                else
                    % There are no general-split constraints.
                    At = zeros([0 q cbSz],'like',inputDataClass);
                    bt = zeros([0 cbSz],'like',inputDataClass);
                    % nrReLUIdx = -ones([0 cbSz],'like',inputDataClass);
                end

                % Refine the input set based on the output specification.
                [li,ui,nrXis] = aux_refineInputSet(nn,options,idxLayer, ...
                    cxi,Gxi,yi,Gyi,A,b,numUnionConst,safeSet,As,bs, ...
                    newNrXi,nrXi,At,bt,nReLU,A_in,b_in);

                % We enclose all unsafe outputs; therefore, a set is
                % verified if it is empty. Identify empty sets. We also
                % flag dummy splits as verified to discard them.
                isDummySplit = repmat(isDummySplit, 1, size(li,2) / size(isDummySplit,2));
                isVerified = any(isnan(li),1) | any(isnan(ui),1) | isDummySplit;
                % Remove the verified sets...
                li(:,isVerified) = [];
                ui(:,isVerified) = [];
                nrXis(:,isVerified) = [];

                % Compute center and radius of refined sets.
                xis = 1/2*(ui + li);
                ris = 1/2*(ui - li);

                % Identify which sets were refined to just being a point.
                isPoint = all(ris == 0,1);
                % Check the specification for the point sets.
                if any(isPoint)
                    [~,~,falsified,x_,y_] = aux_checkPoints(nn,options, ...
                        idxLayer,A,b,safeSet,xis(:,isPoint));
                    if any(falsified)
                        % Validate CE against input-side constraints.
                        if ~isempty(A_in) && ~isempty(x_)
                            if any(A_in * x_ > b_in + 1e-5)
                                falsified = false(size(falsified));
                                x_ = []; y_ = [];
                            end
                        end
                    end
                    if any(falsified)
                        % Found a counterexample.
                        res.str = 'COUNTEREXAMPLE';
                        break;
                    end
                end
                % Remove the point sets from the queue.
                xis(:,isPoint) = [];
                ris(:,isPoint) = [];
                nrXis(:,isPoint) = [];

                % All removed subproblems are verified, i.e., empty set or
                % point sets.
                numVerified = numVerified + sum(isVerified) + sum(isPoint);
                % We have to subtract the number of dummy splits.
                numVerified = numVerified - sum(isDummySplit);
            otherwise
                % Invalid option.
                throw(CORAerror('CORA:wrongFieldValue', ...
                    'options.nn.refinement_method', ...
                    {'naive','zonotack','zonotack-layerwise'}));
        end

        % Add new splits to the queue.
        [xs,rs,nrXs] = aux_push(xis,ris,nrXis,xs,rs,nrXs,options);

        % To save memory, we clear all variables that are no longer used.
        clear(batchVars{:});
        % Clear all batch variables in all layers.
        layers = aux_clearLayerFields(layers);

        % Increment iteration counter.
        iter = iter + 1;
    catch e
        if ismember(e.identifier,{'parallel:gpu:array:pmaxsize', ...
                'parallel:gpu:array:OOM','MATLAB:array:SizeLimitExceeded',...
                'MATLAB:nomem'}) % It is a memory error.
            if nReLU > 0
                % Reduce the number of ReLU constraints.
                if isinf(options.nn.num_relu_constraints)
                    % Reduce to 100 constraints in total.
                    nReLU = 100;
                else
                    % Reduce by factor 10.
                    nReLU = floor(1/10*nReLU);
                end
                fprintf('--- OOM error: reduce ReLU constraints %d...\n',nReLU);
                % Update the options.
                options.nn.num_relu_constraints = nReLU;
            elseif bSz > 1
                % Reduce the batch size.
                bSz = ceil(1/2*bSz);
                fprintf('--- OOM error: half batchSize %d...\n',bSz);
            elseif options.nn.train.num_approx_err > 0
                % Reduce the number of approximation errors.
                if isinf(options.nn.train.num_approx_err)
                    % Reduce the number of approximation errors to 100.
                    options.nn.train.num_approx_err = 1000;
                else
                    % Reduce by factor 10.
                    options.nn.train.num_approx_err = ...
                        floor(1/100*options.nn.train.num_approx_err);
                end
                fprintf('--- OOM error: reduce number of approximation errors %d...\n', ...
                    options.nn.train.num_approx_err);
                % We have to enable the interval center.
                options.nn.interval_center = true;
                % Recompute the indices for the approximation errors.
                q = nn.prepareForZonoBatchEval(x,options,idxLayer);
            elseif numInitGens > 1
                % Reduce the number of input generators by factor 2.
                numInitGens = ceil(1/2*numInitGens);
                fprintf('--- OOM error: reduce number of input generators %d...\n', ...
                    numInitGens);
                % Update the options.
                options.nn.train.num_init_gens = numInitGens;
                % We have to enable the interval center.
                options.nn.interval_center = true;
                % Recompute the indices for the approximation errors.
                q = nn.prepareForZonoBatchEval(x,options,idxLayer);
            else
                % We cannot adapt any parameters.
                fprintf('--- OOM error: no parameters to adapt!\n');
                rethrow(e);
            end
            % Append the current batch items.
            [xs,rs,nrXs] = aux_push(xi,ri,nrXi,xs,rs,nrXs,options);

            % Clear the batch variables.
            clear(batchVars{:});
            % Clear all batch variables in all layers.
            layers = aux_clearLayerFields(layers);
        else
            fprintf('unexpected Error --- neuralNetwork/verify...\n');
            rethrow(e);
        end
    end
end

if size(xs,2) == 0 && ~strcmp(res.str,'COUNTEREXAMPLE')
    % Verified all patches.
    res.str = 'VERIFIED';
    x_ = [];
    y_ = [];
end

% Store time.
res.time = toc(timerVal);
% Store number of verified patches.
res.numVerified = gather(numVerified);

if verbose
    % Print the final stats.
    aux_printIterationStats(table,iter,xs,r,rs,numVerified)
    % Print table footer.
    table.printFooter();
    % Print the result.
    fprintf('--- Result: %s (time: %.3f [s])\n',res.str,res.time);
end

end


% Auxiliary functions -----------------------------------------------------

function [cxi,Gxi,inputDimIds,yi,Gyi,ld_yi,ld_Gyi,ld_Gyi_err,S,sens,grad,a] = ...
    aux_computeBatchArtifacts(nn,options,layers,idxLayer,xi,ri,nrXi,...
    A,b,q,numInitGens,inputGenHeuristic,varargin)
% Compute all per-batch verification artifacts (sensitivity, interval
% gradient, input zonotope, forward propagation, logit difference) for
% the current batch (xi,ri,nrXi).

[computeSens,computeGrads,nNeur,nReLU,nSplits] = ...
    setDefaultValues({false,false,0,0,1},varargin);

% Obtain the current batch size.
[n0,cbSz] = size(xi);

% Compute the sensitivity (used for splitting heuristics, FGSM, etc.).
if computeSens
    [S,~] = nn.calcSensitivity(xi,options,true);
    % Floor sensitivity at 1e-6 to keep splitting heuristic effective.
    sens = reshape(max(max(abs(S),1e-6),[],1),[n0 cbSz]);
else
    S = [];
    sens = [];
end

% Construct the input zonotope.
[cxi,Gxi,inputDimIds] = aux_constructInputZonotope(nn,options,idxLayer, ...
    inputGenHeuristic,xi,ri,A,b,sens,numInitGens,q);

% Apply previously stored neuron splits as layer bounds.
if ~isempty(nrXi)
    layers = aux_setLayerBoundsFromPreviousSplits(layers,nrXi);
end

% Compute and store the gradient that are required for computing the 
% heuristics.
if computeGrads
    grad = aux_updateGradients(nn,options,idxLayer,cxi,Gxi,A,b,false,inputDimIds);
else
    grad = [];
end

% Forward propagation with optional aggregation for split/ReLU
% constraints.
if nNeur > 0 || nReLU > 0
    % Add the function handles to aggregate neuron split and relu
    % constrains.
    options.nn.neuron_aggregation_fun = @(a,layeri,ci,Gi,co,Go) ...
        aux_neuronConstraints(options,nNeur,nSplits,nrXi, ...
        aux_reluConstraints(options,nReLU,nrXi,a,layeri,ci,Gi,co,Go), ...
        layeri,ci,Gi);
    % Do the forward propagation.
    [yi,Gyi,a] = nn.evaluateZonotopeBatch_(cxi,Gxi,options,idxLayer);
else
    % Do the forward propagation.
    [yi,Gyi] = nn.evaluateZonotopeBatch_(cxi,Gxi,options,idxLayer);
    a = struct();
end

% Compute the logit difference.
[ld_yi,ld_Gyi,ld_Gyi_err,~,~,Gyi] = ...
    aux_computeLogitDifference(yi,Gyi,A,options);
end

function aux_printIterationStats(table,iter,xs,r,rs,numVerified)
% Compute the length of the queue.
queueLen = size(xs,2);
% Compute the other stats.
if ~isempty(rs)
    % Compute the average radius.
    avgRad = mean(rs,'all');
    % Compute the ratio of the unknown input volume (assume intervals).
    unknVol = sum(prod(2*rs,1),'all')/sum(prod(2*r),'all')*100;
    % Compute the ratio of f-radii to approximate the volume of unknown
    % input space.
    % unknVol = sum(sqrt(sum(rs.^2,1)),'all')...
    %         /sum(sqrt(sum(r.^2,1)),'all')*100;
else
    % The queue is empty; there are no stats.
    avgRad = 0;
    unknVol = 0;
end
% Print new table row.
table.printContentRow({iter,queueLen,numVerified,avgRad,unknVol});
end

function [xi,ri,nrXi,xs,rs,nrXs,qIdx] = aux_pop(xs,rs,nrXs,bSz,options)
% Obtain the number of elements in the queue.
nQueue = size(xs,2);
% Construct indices to pop.
switch options.nn.verify_dequeue_type
    case 'front'
        % Take the first entries.
        qIdx = 1:min(bSz,nQueue);
    case 'half-half'
        % Half from the front and half from the back.
        qIdx = 1:min(bSz,nQueue);
        offsetIdx = ceil(length(qIdx)/2 + 1):length(qIdx);
        qIdx(offsetIdx) = qIdx(offsetIdx) + nQueue - length(qIdx);
    otherwise
        % Invalid option.
        throw(CORAerror('CORA:wrongFieldValue', ...
            'options.nn.verify_enqueue_type',{'append','prepend'}));
end
% Pop centers.
xi = xs(:,qIdx);
xs(:,qIdx) = [];
% Pop radii.
ri = rs(:,qIdx);
rs(:,qIdx) = [];
% Pop indices for split neurons.
nrXi = nrXs(:,qIdx);
nrXs(:,qIdx) = [];
end

function [xs,rs,nrXs] = aux_push(xi,ri,nrXi,xs,rs,nrXs,options)
% Try to remove irrelevant rows from the stored indices.
nrXi = sort(nrXi,1,'descend'); % Move the NaN to the bottom.
minNumNan = min(sum(isnan(nrXi),1)); % Identify the minimum number of NaN.
nrXi = nrXi((minNumNan+1):end,:); % Trim the minimum number of NaN.

% Retrieve variables from the GPU.
xi = gather(xi);
ri = gather(ri);

% Pad the indices of the split neurons.
if size(nrXi,1) > size(nrXs,1)
    nrXs = [nrXs; NaN([size(nrXi,1)-size(nrXs,1) size(nrXs,2)])];
else
    nrXi = [nrXi; NaN([size(nrXs,1)-size(nrXi,1) size(nrXi,2)])];
end
% Retrieve the indices from the GPU.
nrXi = gather(nrXi);
% Add new splits to the queue.
switch options.nn.verify_enqueue_type
    case 'append'
        % Append the new entries.
        xs = [xs xi];
        rs = [rs ri];
        nrXs = [nrXs nrXi];
    case 'prepend'
        % Prepend the new entries.
        xs = [xi xs];
        rs = [ri rs];
        nrXs = [nrXi nrXs];
    otherwise
        % Invalid option.
        throw(CORAerror('CORA:wrongFieldValue', ...
            'options.nn.verify_enqueue_type',{'append','prepend'}));
end
end

function grad = aux_updateGradients(nn,options,idxLayer,x,rG,A,b, ...
    computeIntervalGradient,inputDimIds)
% Update the gradient of the f-radius store in the layers of the neural
% network. The gradients are used to optimize the approximation slope
% as well as for splitting heuristics.

if computeIntervalGradient
    % Compute the gradient for interval propagation to avoid the full
    % zonotope propagation.
    Yi = nn.evaluate_(interval(x-rG,x+rG),options,idxLayer);

    % We compute the gradient of an interval.
    % Obtain the bounds.
    yli = Yi.inf;
    yui = Yi.sup;

    % Compute the gradient of the f-radius.
    frad = sqrt(sum((yui - yli).^2,1));
    gyli = -yli./(frad + 1e-6);
    gyui = -yui./(frad + 1e-6);
    % Backpropagate the gradients to identify which input dimensions
    % create the largest interval approximation error.
    % For now, we do not store the gradient, because we only use the
    % gradient for the construction of the input zonotope.
    [~,grad] = nn.backpropIntervalBatch(gyli,gyui,options,idxLayer,false);
else
    % Compute the gradient for zonotope propagation.
    [~,Gy] = nn.evaluateZonotopeBatch_(x,rG,options,idxLayer);

    % Obtain the number of output dimensions and batch size.
    [nK,~,bSz] = size(Gy);

    % Compute the gradient of the f-radius.
    frad = sqrt(sum(Gy.^2,[1 2]));
    gGy = Gy./(frad + 1e-6);

    % Compute a dummy center gradient.
    if options.nn.interval_center
        gcy = zeros([nK 2 bSz],'like',Gy);
    else
        gcy = zeros([nK bSz],'like',Gy);
    end

    % Compute gradient of the f-radius of the output set; the gradient
    % is used to split the neuron in the network as well as input
    % dimensions.
    [~,grad] = nn.backpropZonotopeBatch_(gcy,gGy,options,idxLayer,false);
        
    % Obtain the number of input generators.
    [numInitGens,~] = size(inputDimIds);

    % Compute indices for the gradient of the interval
    % norm w.r.t. the different generators.
    dimGenIdx = reshape(sub2ind(size(grad), ...
        inputDimIds, ...
        repmat((1:numInitGens)',1,bSz), ...
        repelem(1:bSz,numInitGens,1)),[numInitGens bSz]);
    % Compute gradient of the norm.
    grad = reshape(grad(dimGenIdx),[numInitGens bSz]);
end
end

function layers = aux_clearLayerFields(layers)
% Remove all stored auxiliary computations in the layers.
for i=1:length(layers)
    % Obtain the i-th layer.
    layeri = layers{i};
    % Clear all fields.
    layeri.backprop.store = struct();
end
end

function [critValPerConstr,critVal,falsified,x_,y_] = ...
    aux_checkPoints(nn,options,idxLayer,A,b,safeSet,xs)
% Compute the output of the adversarial examples.
ys = nn.evaluate_(xs,options,idxLayer);
% Compute the logit difference.
ld_ys = A*ys;
% Check the specification and compute a value indicating how close we
% are to finding an adversarial example (< 0 mean the specification is
% violated).
critValPerConstr = ld_ys - b;
if safeSet
    % safe iff all(A*y <= b) <--> unsafe iff any(A*y > b)
    % Thus, unsafe if any(-A*y < -b).
    falsified = any(ld_ys > b,1);
    critValPerConstr = -critValPerConstr;
    critVal = min(critValPerConstr,[],1);
else
    % unsafe iff all(A*y <= b) <--> safe iff any(A*y > b)
    % Thus, unsafe if all(A*y <= b).
    falsified = all(ld_ys <= b,1);
    critVal = max(critValPerConstr,[],1);
end

if any(falsified)
    % Found a counterexample.
    idNzEntry = find(falsified);
    id = idNzEntry(1);
    x_ = gather(xs(:,id));
    % Gathering weights from gpu. There is are precision error when
    % using single gpuArray.
    nn.castWeights(single(1));
    y_ = nn.evaluate_(x_,options,idxLayer); % yi_(:,id);
else
    % We have not found a counterexample.
    x_ = [];
    y_ = [];
end
end

function [cxi,Gxi,dimIdx] = aux_constructInputZonotope(nn,options,idxLayer, ...
    heuristic,xi,ri,A,b,sens,numInitGens,q)
% Obtain the number of input dimensions and the batch size.
[n0,bSz] = size(xi);

% Initialize the generator matrix.
Gxi = zeros([n0 q bSz],'like',xi);

if numInitGens >= n0
    % We create a generator for each input dimension.
    dimIdx = repmat((1:n0)',1,bSz);
else
    % Compute interval gradient for generator selection.
    if strcmp(heuristic,'zono-norm-gradient')
        grad = aux_updateGradients(nn,options,idxLayer,xi,ri,A,b,true,[]);
    else
        grad = [];
    end

    % Compute the heuristic.
    hi = nnHelper.computeHeuristic(heuristic,xi-ri,xi+ri,ri,sens, ...
        grad,[],[],false,1);

    % Find the input pixels that affect the output the most.
    [~,dimIdx] = sort(hi,'descend');
    % Select the most important input dimensions and add a generator
    % for each of them.
    dimIdx = dimIdx(1:numInitGens,:);
end
% Compute indices for non-zero entries.
gIdx = sub2ind(size(Gxi),dimIdx, ...
    repmat((1:numInitGens)',1,bSz),repelem(1:bSz,numInitGens,1));
% Set non-zero generator entries.
Gxi(gIdx) = ri(sub2ind(size(ri),dimIdx,repelem(1:bSz,numInitGens,1)));
% Sum generators to compute remaining set.
ri_ = (ri - reshape(sum(Gxi,2),[n0 bSz]));

% Construct the center.
if options.nn.interval_center
    % Put remaining set into the interval center.
    cxi = permute(cat(3,xi - ri_,xi + ri_),[1 3 2]);
else
    % The center is just a vector.
    cxi = xi;
end
end

function [l,u] = aux_obtainBoundsFromSplits(neuronIds,bSz,nrSplitIdx,d,c)
% Extract the bounds from constraints. The argument d represents the
% offset of a split constraint; together with the center c we can
% recompute the split-bound.

% Initialize the bounds.
l = -inf([length(neuronIds) bSz],'like',d);
u = inf([length(neuronIds) bSz],'like',d);

% TODO: implement handling of multiple splits per dimension. Use
% aggregation to get the best bounds for each dimension.
if ~isempty(nrSplitIdx) && any(~isnan(nrSplitIdx),'all')
    % Check which constraints contain bounds for the current layer.
    isBnd = (permute(abs(nrSplitIdx),[3 2 1]) ...
        == repmat(neuronIds',1,bSz));
    if ~any(isBnd,'all')
        % There are no bounds from splits.
        return;
    end
    % Check which constraints contain lower or upper bounds.
    isLBnd = (permute(sign(nrSplitIdx),[3 2 1]) == -1);
    isUBnd = (permute(sign(nrSplitIdx),[3 2 1]) == 1);
    % Compute the indices into the bounds.
    ljIdx = any(isBnd & isLBnd,3);
    ujIdx = any(isBnd & isUBnd,3);

    if ~isempty(c) && ~isempty(d)
        % Compute indices into the constraints.
        cstrLjIdx = permute(any(isBnd & isLBnd,1),[3 2 1]);
        cstrUjIdx = permute(any(isBnd & isUBnd,1),[3 2 1]);
        % Extract the bounds based on the splits.
        bndLj = reshape(d(cstrLjIdx),[],1) - c(ljIdx);
        bndUj = reshape(d(cstrUjIdx),[],1) + c(ujIdx);
    else
        % The bounds are 0, because we split a ReLU at 0.
        bndLj = 0;
        bndUj = 0;
    end
    % Update bounds based on the splits.
    if any(ljIdx,'all')
        l(ljIdx) = bndLj;
    end
    if any(ujIdx,'all')
        u(ujIdx) = bndUj;
    end
end
end

function layers = aux_setLayerBoundsFromPreviousSplits(layers,nrXi)
% Use the split indices to set the bounds for previously split ReLU
% neurons.

% Obtain the number of previous split constraints.
[p,bSz] = size(nrXi);
% For now we only support storing ReLU splits at 0.
% TODO: remember arbitrary splits.
% Therefore, create dummy offsets.
dummyd = zeros([p bSz],'like',nrXi);
% Store computed bounds in the layers for tighter
% approximations.
for i=1:length(layers)
    % Obtain the i-th layer.
    layeri = layers{i};
    if isa(layeri,'nnActivationLayer')
        % Obtain the indices of the neurons of the current
        % layer.
        neuronIds = layeri.neuronIds;
        % Create dummy centers.
        dummyc = zeros([length(neuronIds) bSz],'like',nrXi);
        % Compute bounds from previous splits.
        [li,ui] = aux_obtainBoundsFromSplits(neuronIds,bSz, ...
            nrXi,dummyd,dummyc);
        if any(~isinf(li) | ~isinf(ui),'all')
            % Store the computed bounds in the layers.
            layeri.backprop.store.l = li;
            layeri.backprop.store.u = ui;
        end
    end
end
end

function [l,u,cl,cu] = aux_computeBoundsOfZonotope(c,G,options)
% Compute the bounds of a batch of zonotopes. We return the bounds of
% the refined set, the unrefined (only duplicated set), as well as the
% refined set.

% Obtain number of hidden neurons.
[nk,~,bSz] = size(G);

% Compute the radius of the zonotope.
r = reshape(sum(abs(G),2),[nk bSz]);
if options.nn.interval_center
    % Compute center and center radius.
    cl = reshape(c(:,1,:),[nk bSz]);
    cu = reshape(c(:,2,:),[nk bSz]);
else
    % The radius is zero.
    cl = c;
    cu = c;
end
% Compute the bounds.
l = cl - r;
u = cu + r;
end

function a = aux_updateLayerBounds(options,numUnionConst,C,d,a,layeri,ci,Gi)
% (Not nice!) Abuse neuron aggregation to set layer bounds from a
% restricted factor space.

if isa(layeri,'nnActivationLayer')
    % We only need to set the bounds for activation layers.

    % Obtain the number of dimensions, number of generators, and batch
    % size.
    [ni,~,bSz] = size(Gi);

    if options.nn.interval_center
        % Compute center and center radius.
        cr = reshape(1/2*(ci(:,2,:) - ci(:,1,:)),[ni bSz]);
        ci = reshape(1/2*(ci(:,2,:) + ci(:,1,:)),[ni bSz]);
    else
        % The radius is zero.
        cr = 0;
    end

    % Construct a struct for the output set.
    uXi = struct('c',ci,'dr',cr,'G',Gi);
    % Apply the constraints.
    uXi.A = C;
    uXi.b = d;

    % Compute the bounds of the unsafe inputs (hypercube).
    [li,ui,~,~] = conZonotope.approximateBoundsWithGPU( ...
        uXi,numUnionConst,options);

    % Set the bounds for the next propagation.
    layeri.backprop.store.l_ = li;
    layeri.backprop.store.u_ = ui;
end
end

function [ld_yi,ld_Gyi,ld_Gyi_err,yic,yid,Gyi] = ...
    aux_computeLogitDifference(yi,Gyi,A,options)
% Obtain number of output dimensions and batch size.
[nK,~,bSz] = size(Gyi);

if options.nn.interval_center
    % Compute the center and the radius of the center-interval.
    yic = reshape(1/2*(yi(:,2,:) + yi(:,1,:)),[nK bSz]);
    % Compute approximation error.
    yid = 1/2*(yi(:,2,:) - yi(:,1,:));
else
    % The center is just a vector.
    yic = yi;
    % There are no approximation errors stored in the center.
    yid = zeros([nK 1 bSz],'like',yi);
end

% Compute the logit difference of the input generators.
ld_yi = A*yic;
ld_Gyi = pagemtimes(A,Gyi);
% Compute logit difference of the approximation errors.
ld_Gyi_err = sum(abs(A.*permute(yid,[2 1 3])),2);
end

function [c,G] = aux_matchBatchSize(c,G,bSz)
% Identify if c is an interval center (iff options.nn.interval_center).
intervalCenter = size(c,3) > 1 || (size(c,2) == 2 && size(G,3) == 1);
% Replicate a zonotope batch for splitting.
if bSz ~= size(G,3) % iff newSplits > 1
    newSplits = bSz/size(G,3);
    if intervalCenter
        c = repelem(c,1,1,newSplits);
    else
        c = repelem(c,1,newSplits);
    end
    G = repelem(G,1,1,newSplits);
end
end

% Set Refinement ----------------------------------------------------------

function [l,u,nrXis] = aux_refineInputSet(nn,options,idxLayer, ...
    x,Gx,y,Gy,A,b,numUnionConst,safeSet,As,bs,newNrXs,prevNrXs,At,bt,nReLU, ...
    A_in,b_in)

% Specify the maximum number of refinement iterations.
maxRefIter = options.nn.refinement_max_iter;

% Handle optional input-side linear constraints (A_in*x <= b_in).
if nargin < 20, A_in = []; b_in = []; end

% Obtain number of input dimensions, generators, and batchsize.
n0 = size(x,1);
[~,~,bSz] = size(Gy);

% Convert and join the general- & input-split constraints.
[C,d,newSplits,nrSplitIdx] = aux_convertSplitConstraints(As,bs,newNrXs);

if ~isempty(At) && ~isempty(bt)
    % Replicate the ReLU constraints to match the batch size.
    [bt,At] = aux_matchBatchSize(bt,At,bSz*newSplits);
end

% Replicate set for split constraints.
[x,Gx] = aux_matchBatchSize(x,Gx,bSz*newSplits);
[y,Gy] = aux_matchBatchSize(y,Gy,bSz*newSplits);
% Update the batch size.
bSz = bSz*newSplits;

if ~isempty(nrSplitIdx) && size(bs,2) == 1
    % Compute the indices of the split neurons.
    nrXis = nrSplitIdx;
    % Remove all indices that do not correspond to the split of a neuron.
    nrXis(all(isnan(nrXis),2),:) = [];
    nrXis = reshape(nrXis,[],bSz);
    % Duplicate the indices for previously split neurons.
    prevNrXs = repelem(prevNrXs,1,newSplits);
    % Combine the newly split neurons with the previously split ones.
    nrXis = [prevNrXs; nrXis];
else
    % Currently, we can only remember splits into two pieces.
    nrXis = zeros([0 bSz]);
end

% Enumerate the layers of the neural networks.
layers = nn.enumerateLayers();
if ~isempty(nrXis)
    % Set layer bounds based on previously split ReLU neurons.
    layers = aux_setLayerBoundsFromPreviousSplits(layers,nrXis);
end

% Initialize scale and offset of the generators.
bc = 0;
br = 1;

% Obtain the input generators indices.
qiIds = nn.layers{1}.genIds;

% Initialize the counter for number of refinement iterations.
refIter = 1;

% Keep track of empty sets.
isEmpty = zeros([1 bSz],'logical');

% Iteratively refine the input set to enclose only the unsafe outputs.
while refIter < maxRefIter

    % Construct the unsafe input set.
    uXi = aux_constructUnsafeInputSet(options,x,Gx,qiIds,y,Gy,A,b, ...
        safeSet,numUnionConst);

    if ~isempty(C)
        % Scale and offset constraints with current hypercube.
        [d,C] = aux_scaleAndOffsetZonotope(d,C,-bc,br);

        % Append split constraints.
        uXi.A = [uXi.A; C];
        uXi.b = [uXi.b; d];

        if ~isempty(At) && ~isempty(bt)
            % Scale and offset constraints with current hypercube.
            [bt,At] = aux_scaleAndOffsetZonotope(bt,At,-bc,br);
            % Append the ReLU tightening constraints.
            uXi.A = [uXi.A; At];
            uXi.b = [uXi.b; bt];
        end
    end

    % Inject input-side polytope constraints: A_in*(xc+Gx*beta) <= b_in.
    if ~isempty(A_in)
        qGx  = size(Gx,2);
        bSzG = size(Gx,3);
        if options.nn.interval_center && ndims(x) > 2
            xc_ = reshape(1/2*(x(:,2,:) + x(:,1,:)),[n0 bSzG]);
        else
            xc_ = x;
        end
        A_in_s = cast(A_in,'like',Gx);
        b_in_s = cast(b_in,'like',Gx);
        pIn  = size(A_in_s,1);
        % (p_in × n0) * (n0 × q×bSzG) -> (p_in × q × bSzG)
        C_in = reshape(A_in_s * reshape(Gx,[n0, qGx*bSzG]),[pIn, qGx, bSzG]);
        % (p_in × 1) broadcast against (p_in × bSzG) -> (p_in × bSzG)
        d_in = b_in_s - A_in_s * xc_;
        uXi.A = [uXi.A; C_in];
        uXi.b = [uXi.b; d_in];
    end

    % Compute the bounds of the unsafe inputs (hypercube).
    [li_,ui_,bli,bui] = conZonotope.approximateBoundsWithGPU( ...
        uXi,numUnionConst,options);
    % Convert bounds of the hypercube to center and radius.
    br = 1/2*(bui - bli);
    bc = 1/2*(bui + bli);

    % Update empty sets.
    isEmpty = isEmpty | any(isnan(li_),1) | any(isnan(ui_),1);

    % Scale and offset the input set.
    [x,Gx] = aux_scaleAndOffsetZonotope(x,Gx,bc,br);

    % Increment the refinement iteration counter.
    refIter = refIter + 1;

    if refIter < maxRefIter
        % Compute a new output enclosure.
        [y,Gy] = nn.evaluateZonotopeBatch_(x,Gx,options,idxLayer);
    end
end

% Compute bounds of the refined input set.
[l,u] = aux_computeBoundsOfZonotope(x,Gx,options);

% Clear all variables in all layers.
layers = aux_clearLayerFields(layers);

% Update empty sets.
isEmpty = isEmpty | any(l > u,1);

% Update bounds to represent empty sets.
l(:,isEmpty) = NaN;
u(:,isEmpty) = NaN;
end

function uX = aux_constructUnsafeInputSet(options,x,Gx,qIds, ...
    y,Gy,A,b,safeSet,numUnionConst)

% Obtain the number of input dimensions, input generators, and batch size.
[n0,qx,bSz] = size(Gx);
% Obtain the number of output generators.
[~,qy,~] = size(Gy);
% Compute the generator indices for the approximaion errors.
% eIds = setdiff(1:qy,qIds);
% Obtain the number of output constraints.
% [p,~,~] = size(A);

if options.nn.interval_center
    % Compute center and center radius.
    xc = reshape(1/2*(x(:,2,:) + x(:,1,:)),[n0 bSz]);
    xr = reshape(1/2*(x(:,2,:) - x(:,1,:)),[n0 bSz]);
else
    % The radius is zero.
    xc = x;
    xr = 0;
end

% Compute the output constraints.
[ld_yi,ld_Gyi,ld_Gyi_err,~,~,~] = ...
    aux_computeLogitDifference(y,Gy,A,options);
% Compute output constraints.
if safeSet
    % safe iff all(A*y <= b) ...
    % <--> unsafe iff any(A*y > b) <--> unsafe iff any(-A*y < -b)
    % Thus, unsafe if any(-A*Gy*\beta < -b + A*y)
    A_ = -ld_Gyi;
    b_ = ld_yi - b;
    % Invert the sign for the union constraints.
    A_((numUnionConst+1):end,:,:) = -A_((numUnionConst+1):end,:,:);
    b_((numUnionConst+1):end,:) = -b_((numUnionConst+1):end,:);
else
    % unsafe iff all(A*y <= b)
    % Thus, unsafe if all(A*Gy*\beta <= b - A*y)
    A_ = ld_Gyi;
    b_ = b - ld_yi;
end

% Construct a struct for the output set.
uX = struct('c',xc,'dr',xr,'G',Gx);
% Apply the output constraints to the input set of the i-th layer.
uX.A = A_; % zeros([p qx bSz],'like',Gx);
% uX.A(:,qIds,:) = A_(:,qIds,:);
% Offset by refinement errors.
uX.b = b_ + ld_Gyi_err(:,:);  ...
    % + reshape(sum(abs(A_(:,eIds,:)),2),[p bSz]);
end

function [c,G] = aux_scaleAndOffsetZonotope(c,G,bc,br)
% Obtain indices of generator.
qiIds = 1:min(size(G,2),size(bc,1));
% Extract the relevant entries.
G_ = G(:,qiIds,:);
bc_ = permute(bc(qiIds,:),[1 3 2]);
br_ = permute(br(qiIds,:),[3 1 2]);
% Scale and offset the zonotope to a new hypercube with center bic and
% radius bir.
offset = pagemtimes(G_,bc_);
% Offset the center.
if ndims(c) > 2 % iff options.nn.interval_center
    c = c + offset;
else
    c = c + offset(:,:);
end
% Scale the generators.
G(:,qiIds,:) = G(:,qiIds,:).*br_;
end

function [A,b,newSplits,nrSplitIdx] = aux_convertSplitConstraints( ...
    As,bs,nrXis)
% Consider all combinations between the given constraints, A*x <= b.
if ~isempty(As)
    % Obtain the number of split-constraints.
    [ps,q,bSz] = size(As);
    % Obtain the number of pieces.
    [~,pcs,~] = size(bs);
    % Compute number of new splits.
    newSplits = (pcs+1)^ps;

    % We flip the signs of the constraints to realize splitting; -1
    % represents a lower bound, i.e., -A*x > -b, whereas 1 represents
    % an upper bound, i.e., A*x <= b.
    constrSign = [-1; 1];

    % Duplicate each halfspace for a lower and an upper bound.
    As_ = repelem(As,2,1,1);
    % Duplicate offsets for lower and upper bound.
    bs_ = repelem(bs,1,2,1);
    % Duplicate the indices for the split neurons.
    nrSplitIdx = repelem(nrXis,2,1);

    % Scale the constraints; -1 for upper bound and 1 for lower bound.
    As_ = repmat(constrSign,ps,1).*As_;
    % Duplicate the constraint for the new splits.
    A_ = permute(repelem(As_,1,1,1,newSplits),[2 1 4 3]);

    % Mark unused bounds by NaN.
    bs_ = cat(2,NaN(ps,1,bSz),bs_,cat(2,NaN(ps,1,bSz)));
    % Scale the offsets; -1 for upper bound and 1 for lower bound.
    bs_ = repmat(constrSign',1,pcs+1).*bs_;
    % Reshape and combine the lower and upper bounds.
    bs_ = reshape(permute(reshape(permute(bs_,[4 2 1 3]), ...
        [2 pcs+1 ps bSz]),[1 3 2 4]),[2*ps pcs+1 bSz]);
    % Extend the offsets.
    b_ = [bs_ zeros([2*ps newSplits - (pcs+1) bSz],'like',bs_)];
    % Compute all combinations of the splits.
    idx = pcs+1;
    for i=1:(ps-1)
        % Increase the index.
        idx_ = idx*(pcs+1);
        % Repeat the current combined splits.
        b_(1:2*i,1:idx_,:) = repmat(b_(1:2*i,1:idx,:),1,pcs+1,1);
        % Repeat the elements of the next split and append them.
        b_(2*i + (1:2),1:idx_,:) = ...
            repelem(b_(2*i + (1:2),1:(pcs+1),:),1,(pcs+1)^i,1);
        % Update the index of the combined splits.
        idx = idx_;
    end

    % Compute the neuron indices of each constraint. Then we can track
    % from which neuron the constraint stems; to extract the exact
    % bounds of the set.
    % Scale the indices; -1 for upper bound and 1 for lower bound.
    nrSplitIdx = repmat(constrSign,ps,1).*nrSplitIdx;
    % Duplicate the indices for the new splits.
    nrSplitIdx = repelem(nrSplitIdx,1,newSplits);

    % Find all unused constraints.
    nanIdx = isnan(b_);
    % Set all not needed constraints to zero.
    A_(:,nanIdx) = 0;
    b_(nanIdx) = 0;
    % Mark all unused constraints.
    nrSplitIdx(nanIdx) = NaN;

    % Reshape the constraint matrix and offset.
    A = reshape(permute(A_,[2 1 3 4]),[2*ps q newSplits*bSz]);
    b = reshape(b_,[2*ps newSplits*bSz]);
else
    % There are no additional constraints.
    newSplits = 1;
    A = zeros([0 size(As,[2 3])],'like',As);
    b = zeros([0 size(bs,2)],'like',bs);
    nrSplitIdx = zeros([0 size(nrXis,2)],'like',nrXis);
end
end

% Constraints & Splitting -------------------------------------------------

function [xis,ris,dimId] = aux_split(xi,ri,hi,nSplits)
% Split one input dimension into nSplits pieces.
[n,bSz] = size(xi);
% Split each input in the batch into nSplits parts.
% 1. Find the input dimension with the largest heuristic.
[~,sortDims] = sort(hi,1,'descend');
dimIds = sortDims(1,:);
% Construct indices to use sub2ind to compute the offsets.
splitsIdx = repmat(1:nSplits,1,bSz);
bSzIdx = repelem((1:bSz)',nSplits);

dimId = dimIds(1,:);
linIdx = sub2ind([n bSz nSplits], ...
    repelem(dimId,nSplits),bSzIdx(:)',splitsIdx(:)');
% 2. Split the selected dimension.
xi_ = xi;
ri_ = ri;
% Shift to the lower bound.
dimIdx = sub2ind([n bSz],dimId,1:bSz);
xi_(dimIdx) = xi_(dimIdx) - ri(dimIdx);
% Reduce radius.
ri_(dimIdx) = ri_(dimIdx)/nSplits;

xis = repmat(xi_,1,1,nSplits);
ris = repmat(ri_,1,1,nSplits);
% Offset the center.
xis(linIdx(:)) = xis(linIdx(:)) + (2*splitsIdx(:) - 1).*ris(linIdx(:));

% Flatten.
xis = xis(:,:);
ris = ris(:,:);
% Replicate the dimension index.
dimId = repmat(dimId,1,nSplits);
end

function [Ai,bi,dimIds,hi] = aux_dimSplitConstraints(hi,nSplits,nDims)
% Construct dimension split constraints that splits #nDims dimensions
% into #nSplits pieces for subsequent refinement.

% Obtain the number of dimensions and batch size.
[n,bSz] = size(hi);
nDims = min(nDims,n);

% Split each input in the batch into nSplits parts.
% 1. Find the input dimension with the largest heuristic.
[hi,sortDims] = sort(hi,1,'descend');
dimIds = sortDims(1:nDims,:);
hi = hi(1:nDims,:);

% Compute dimension indices.
dimIdx = sub2ind([nDims n bSz],repelem((1:nDims)',1,bSz), ...
    dimIds,repelem(1:bSz,nDims,1));

% 2. Construct the constraints.
Ai = zeros([nDims n bSz],'like',hi);
Ai(dimIdx) = 1;
bi = repelem(-1 + (1:(nSplits-1)).*(2/nSplits),nDims,1,bSz); % Specify offsets.
end

function [er,sens,grad,neuronIds] = ...
    aux_extractLayerFieldsForHeuristicComputation(layeri,ni,bSz)
% Obtain the indices of the neurons of the current layer.
neuronIds = layeri.neuronIds;

% Obtain the approximation errors.
if isfield(layeri.backprop.store,'el') && ...
        isfield(layeri.backprop.store,'eu')
    el = layeri.backprop.store.el;
    eu = layeri.backprop.store.eu;
else
    % We have not stored the approximation errors.
    el = [];
    eu = [];
end
% Compute the radius of approximation errors.
er = 1/2*(eu - el);
er = aux_sliceToBatch(er,bSz);

% Obtain the sensitivity for heuristic.
if ~isempty(layeri.sensitivity)
    Si_ = max(abs(layeri.sensitivity),1e-6);
    % When cbSz==1 MATLAB stores sensitivity as 2-D [nK,ni]; size(.,3)==1.
    % Slice or replicate to match the current bSz.
    if size(Si_,3) > bSz
        Si_ = Si_(:,:,1:bSz);
    elseif size(Si_,3) < bSz
        Si_ = repmat(Si_,[1 1 ceil(bSz/size(Si_,3))]);
        Si_ = Si_(:,:,1:bSz);
    end
    % Max over output dimensions and flatten to [ni, bSz].
    sens = reshape(max(Si_,[],1),ni,bSz);
else
    sens = [];
end

% Obtain the gradient of the zonotope norm.
if isfield(layeri.backprop.store,'approx_error_gradients')
    % Obtain the stored gradient.
    grad = layeri.backprop.store.approx_error_gradients;
    grad = aux_sliceToBatch(grad,bSz);
else
    % There is no stored gradient.
    grad = [];
end
end

function M = aux_sliceToBatch(M,bSz)
% Slice the last dimension of M to the first bSz entries. Input
% splitting can reduce the current bSz below the batch size at which
% these tensors were originally stored.
if isempty(M)
    return;
end
nd = ndims(M);
curBSz = size(M,nd);
if curBSz == bSz
    return;
end
if curBSz < bSz
    reps = ones(1,nd);
    reps(nd) = ceil(bSz/curBSz);
    M = repmat(M,reps);
end
if nd == 2
    M = M(:,1:bSz);
else
    idx = repmat({':'},1,nd);
    idx{nd} = 1:bSz;
    M = M(idx{:});
end
end

function a = aux_neuronConstraints(options,nNeur,nSplits,prevNrXs,a,layeri,ci,Gi)
% Function handle for neuron aggregation to construct the neuron split
% constraints.

if nNeur == 0
    return;
end

% Obtain the number of dimensions, number of generators, and batch
% size.
[ni,qi,bSz] = size(Gi);

% Initialize the aggregation result.
if ~isfield(a,'As') || ~isfield(a,'bs') ...
        || ~isfield(a,'h') || ~isfield(a,'nrSplitIdx')
    % Initialize constraints.
    a.As = zeros([nNeur qi bSz],'like',ci);
    a.bs = zeros([nNeur nSplits-1 bSz],'like',ci);
    % Initial heuristics.
    a.h = -ones([nNeur bSz],'like',ci);
    % Initialize indices of neuron split.
    a.nrSplitIdx = Inf([nNeur bSz],'like',ci);
end

if isa(layeri,'nnActivationLayer')
    % We only generate splits for the input set of activation layers.

    % Compute the bounds of the input set.
    [li,ui] = aux_computeBoundsOfZonotope(ci,Gi,options);
    % Extract the relevent layer fields for the heuristics computation.
    [er,sens,grad,neuronIds] = ...
        aux_extractLayerFieldsForHeuristicComputation(layeri,ni,bSz);
    % Extract the heuristic for splitting neurons.
    heuristic = options.nn.neuron_split_heuristic;
    % Compute the heuristic.
    hi = nnHelper.computeHeuristic(heuristic,li,ui,er,sens,grad,prevNrXs,neuronIds);

    % Create new constraints, i.e., Asi*beta<=bsi.
    Ai = Gi;
    switch options.nn.neuron_split_position
        case 'zero'
            % Split into #nSplits pieces around 0.
            nSplits_ = floor((nSplits-1)/2);
            splitEnum = 1/(nSplits_+1).*(1:floor((nSplits-1))/2)';
            bil = flip(splitEnum).*permute(li,[3 1 2]);
            biu = splitEnum.*permute(ui,[3 1 2]);
            if mod(nSplits,2) == 0
                % Include the center in the lower bounds.
                bil = [bil; zeros([1 ni bSz],'like',ci)];
            end
            % Combine the bounds.
            bi = [bil; biu];
            % Compute the center (ci could be an interval center).
            ci_ = 1/2*(ui + li);
            % Subtract the center.
            bi = reshape(bi - permute(ci_,[3 1 2]),[nSplits-1 ni bSz]);
        case 'middle'
            % Split into #nSplits pieces around the middle.
            splitEnum = linspace(-1,1,nSplits+1)';
            splitEnum = splitEnum(2:end-1);
            % Compute the radius.
            ri = 1/2*(ui - li);
            % Scale the radius.
            bi = splitEnum.*permute(ri,[3 1 2]);
        otherwise
            % Invalid option.
            throw(CORAerror('CORA:wrongFieldValue', ...
                'options.nn.neuron_split_position',{'zero','middle'}));
    end
    % Transpose the constraint offsets.
    bi = permute(bi,[2 1 3]);

    % Select the top constraints for splitting.
    [a.As,a.bs,a.h,idx] = aux_extractTopKConstraints( ...
        nNeur,a.As,a.bs,a.h,Ai,bi,hi);

    % Update the splitting indices.
    nrSplitIdx = [a.nrSplitIdx; repmat(neuronIds,bSz,1)'];
    a.nrSplitIdx = reshape(nrSplitIdx(idx),[nNeur bSz]);
end
end

function a = aux_reluConstraints(options,nReLU,prevNrXs,a,layeri,ci,Gi,co,Go)
% Function handle for neuron aggregation to construct ReLU constraints:
% (i)   ReLU(x) >= 0,
% (ii)  ReLU(x) >= x, and
% (iii) ReLU(x) <= u/(u-l)*(x - l).

if nReLU == 0
    return;
end

% Obtain the number of dimensions, number of generators, and batch
% size.
[ni,qi,bSz] = size(Gi);

% Initialize the aggregation result.
if ~isfield(a,'At1') || ~isfield(a,'bt1') || ...
        ~isfield(a,'At2') || ~isfield(a,'bt2') || ...
        ~isfield(a,'ht') || ~isfield(a,'nrReLUIdx')
    % Initialize constraints.
    a.At1 = zeros([nReLU qi bSz],'like',ci);
    a.bt1 = zeros([nReLU 1 bSz],'like',ci);
    a.At2 = zeros([nReLU qi bSz],'like',ci);
    a.bt2 = zeros([nReLU 1 bSz],'like',ci);
    % Initial heuristics.
    a.ht = -ones([nReLU bSz],'like',ci);
    % Initialize indices of neurons.
    a.nrReLUIdx = Inf([nReLU bSz],'like',ci);
end

if isa(layeri,'nnActivationLayer')
    % We only generate splits for the input set of activation layers.

    % Compute the bounds of the input set.
    [li,ui,cil,~] = aux_computeBoundsOfZonotope(ci,Gi,options);
    % Extract the relevent layer fields for the heuristics computation.
    [er,sens,grad,neuronIds] = ...
        aux_extractLayerFieldsForHeuristicComputation(layeri,ni,bSz);
    % Extract the heuristic for ReLU constraints.
    heuristic = options.nn.relu_constraint_heuristic;
    % Compute the heuristic.
    hi = nnHelper.computeHeuristic(heuristic,li,ui,er,sens,grad,prevNrXs,neuronIds);

    % Compute the upper center bound of the output set.
    [~,~,~,co_u] = aux_computeBoundsOfZonotope(co,Go,options);

    % Create new constraints.

    % (i) ReLU(x) >= 0
    % y_j = co_j + Go_j*\beta >= 0  <-->  -Go_j*\beta <= co_j
    Ai1 = -Go;
    bi1 = permute(co_u,[1 3 2]);

    % (ii) ReLU(x) >= x
    % y_j >= x_j  <-->  (Gi_j - Go_j)*\beta <= co_j - ci_j
    Ai2 = Gi - Go;
    bi2 = permute(co_u - cil,[1 3 2]);

    % Select the top constraints.
    [a.At1,a.bt1,~,~] = aux_extractTopKConstraints( ...
        nReLU,a.At1,a.bt1,a.ht,Ai1,bi1,hi);
    [a.At2,a.bt2,a.ht,idx] = aux_extractTopKConstraints( ...
        nReLU,a.At2,a.bt2,a.ht,Ai2,bi2,hi);

    % Update the splitting indices.
    nrReLUIdx = [a.nrReLUIdx; repmat(neuronIds,bSz,1)'];
    a.nrReLUIdx = reshape(nrReLUIdx(idx),[nReLU bSz]);
end
end

function [At,bt] = aux_extractAndReduceReluConstraints(a)
% Create tightening constraints for the ReLU enclosures.

% Extract the computed constraints.
At = [a.At1; a.At2];
bt = [a.bt1; a.bt2];
% Compute the number of different ReLU constraints.
% numTypesOfConstr = size(At,1)/size(a.nrReLUIdx,1);
% nrReLUIdx = repmat(a.nrReLUIdx,numTypesOfConstr,1);

% Obtain the batch size.
[p,~,bSz] = size(At);
% Reshape the offset.
bt = reshape(bt,[p bSz]);

% Sort the constraints and to remove invalid constraints.
[~,sortIds] = sort(isnan(bt),1,'descend');
% Find the minimal number of invalid constraints across the batch.
minNumInvalidConstr = min(sum(any(isnan(At),2),1),[],3);
% Compute linear indices to remove the constraints.
sortIdx = sub2ind(size(At,[1 3]), ...
    sortIds,repmat(1:bSz,[p 1]));
% Reorder the constraints s.t. the invalid constraints
% are on top.
At = permute(At,[2 1 3]);
At = reshape(At(:,sortIdx),size(At));
bt = reshape(bt(sortIdx),size(bt));
% nrReLUIdx = reshape(nrReLUIdx(sortIdx),size(nrReLUIdx));
% Remove the invalid constraints.
At(:,1:minNumInvalidConstr,:) = [];
bt(1:minNumInvalidConstr,:) = [];
% nrReLUIdx(1:minNumInvalidConstraints,:) = [];

% Set all remaining invalid constraints to 0.
At(isnan(At)) = 0;
bt(isnan(bt)) = 0;

% Transpose constraint matrix.
At = permute(At,[2 1 3]);
end

function [A,b,h,idx] = aux_extractTopKConstraints(k,A1,b1,h1,A2,b2,h2)
% Given two set of constraints with heuristic values, we select the top
% k.

% Concatenate the constraints.
A = [A1; A2];
b = [b1; b2];

% Obtain the number of dimensions and batch size.
[~,q,bSz] = size(A);
% Obtain the number of splits.
[~,s,~] = size(b);

% Sort the heuristic to identify the best constraints.
[h,ids] = sort([h1; h2],1,'descend');
% Only keep the constraints for the top neurons.
h = h(1:k,:);
% Obtain the indices for the relevant constraints.
idx = sub2ind(size(A,[1 3]),ids(1:k,:),repmat(1:bSz,k,1));

% Use indexing to extract the selected constraints.
A = permute(A,[2 1 3]);
A = permute(reshape(A(:,idx),[q k bSz]),[2 1 3]);
b = permute(b,[2 1 3]);
b = permute(reshape(b(:,idx),[s k bSz]),[2 1 3]);
end

% ------------------------------ END OF CODE ------------------------------
