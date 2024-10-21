function [res, x_, y_] = verify(nn, x, r, A, b, safeSet, options, timeout, verbose)
% verify - automated verification for specification on neural networks.
%
% Syntax:
%    [res, z] = nn.verify(x, r, A, b, options)
%
% Inputs:
%    nn - object of class neuralNetwork
%    x, r - center and radius of the initial set
%    A, b - specification, prove A*y + b <= 0
%    safeSet - bool, safe-set or unsafe-set
%    options - evaluation options
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

% Authors:       Niklas Kochdumper, Tobias Ladner, Lukas Koller
% Written:       23-November-2021
% Last update:   30-November-2022 (TL, removed neuralNetworkOld, adaptive)
%                25-July-2023 (TL, input parsing, improvements)
%                23-November-2023 (TL, verbose, bug fix)
%                14-June-2024 (LK, rewritten with efficient splitting)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

nSplits = 5; % 5;
nDims = 1; % 5;

totalNumSplits = 0;
verifiedPatches = 0;

% Extract parameters.
bs = options.nn.train.mini_batch_size;

% To speed up computations and reduce gpu memory, we only use single 
% precision.
inputDataClass = single(1);
% Check if a gpu is used during training.
useGpu = options.nn.train.use_gpu;
if useGpu
    % Training data is also moved to gpu.
    inputDataClass = gpuArray(inputDataClass);
end
% (potentially) move weights of the network to gpu
nn.castWeights(inputDataClass);

% Specify indices of layers for propagation.
idxLayer = 1:length(nn.layers);

% In each layer, store ids of active generators and identity matrices 
% for fast adding of approximation errors.
numGen = nn.prepareForZonoBatchEval(x,options,idxLayer);
% Allocate generators for initial perturbance set.
idMat = cast([eye(size(x,1)) zeros(size(x,1),numGen - size(x,1))], ...
    'like',inputDataClass);
batchG = cast(repmat(idMat,1,1,bs),'like',inputDataClass);

% Initialize queue.
xs = x;
rs = r;
% Obtain number of input dimensions.
n0 = size(x,1);

res = [];

tic

% Main splitting loop.
while size(xs,2) > 0

    time = toc;
    if time > timeout
        res = 'UNKNOWN';
        x_ = [];
        y_ = [];
        break;
    end

    if verbose
        fprintf('Queue / Verified / Total: %07d / %07d / %07d [Avg. radius: %.5f]\n', ...
            size(xs,2),verifiedPatches,totalNumSplits,mean(rs,'all'));
    end

    % Pop next batch from the queue.
    [xi,ri,xs,rs] = aux_pop(xs,rs,bs);
    % Move the batch to the GPU.
    xi = cast(xi,'like',inputDataClass);
    ri = cast(ri,'like',inputDataClass);

    % Falsification -------------------------------------------------------
    % Try to falsification with a FGSM attack.
    % 1. Compute the sensitivity.
    [S,~] = nn.calcSensitivity(xi,options,false);
    S = max(S,1e-3);
    % sens = permute(sum(pagemtimes(A,S)),[2 1 3]);
    sens = permute(sum(abs(S)),[2 1 3]);
    sens = sens(:,:);
    % 2. Compute adversarial attacks. We want to maximze A*yi + b; 
    % therefore, ...
    zi = xi + ri.*sign(sens);
    % 3. Check adversarial examples.
    yi = nn.evaluate_(zi,options,idxLayer);
    if safeSet
        checkSpecs = any(A*yi + b >= 0,1);
    else
        checkSpecs = all(A*yi + b <= 0,1);
    end
    if any(checkSpecs)
        % Found a counterexample.
        res = 'COUNTEREXAMPLE';
        idNzEntry = find(checkSpecs);
        id = idNzEntry(1); % ceil(idNzEntry/size(A,1));
        x_ = zi(:,id);
        % Gathering weights from gpu. There is are precision error when 
        % using single gpuArray.
        nn.castWeights(single(1));
        y_ = nn.evaluate_(gather(x_),options,idxLayer); % yi(:,id);
        break;
    end

    % Verification --------------------------------------------------------
    % 1. Use batch-evaluation.
    if ~options.nn.interval_center
        cxi = xi;
    else
        cxi = permute(repmat(xi,1,1,2),[1 3 2]);
    end
    Gxi = permute(ri,[1 3 2]).*batchG(:,:,1:size(ri,2));
    [yi,Gyi] = nn.evaluateZonotopeBatch_(cxi,Gxi,options,idxLayer);
    % 2. Compute logit-difference.
    if ~options.nn.interval_center
        dyi = A*yi + b;
        dri = sum(abs(pagemtimes(A,Gyi)),2);
    else
        % Compute the center and the radius of the center-interval.
        yic = 1/2*(yi(:,2,:) + yi(:,1,:));
        yid = 1/2*(yi(:,2,:) - yi(:,1,:));
        % Comute the logit difference.
        dyi = A*yic(:,:) + b;
        dri = sum(abs(pagemtimes(A,Gyi)),2) ...
            + sum(abs(A.*pagetranspose(yid)),2);
    end
    % 3. Check specification.
    if safeSet
        checkSpecs = any(dyi(:,:) + dri(:,:) > 0,1);
    else
        checkSpecs = all(dyi(:,:) - dri(:,:) < 0,1);
    end
    unknown = checkSpecs;
    xi = gather(xi);
    ri = gather(ri);
    sens = gather(sens);
    % 3. Create new splits.
    [xis,ris] = aux_split(xi(:,unknown),ri(:,unknown),sens(:,unknown), ...
        nSplits,nDims);
    % Add new splits to the queue.
    xs = [xis xs];
    rs = [ris rs];

    totalNumSplits = totalNumSplits + size(xis,2);
    verifiedPatches = verifiedPatches + size(xi,2) - sum(unknown,'all');

    % To save memory, we clear all variables that are no longer used.
    batchVars = {'xi','ri','xGi','yi','Gyi','dyi','dri'};
    clear(batchVars{:});
end

% Verified.
if isempty(res)
    res = 'VERIFIED';
    x_ = [];
    y_ = [];
end

end


% Auxiliary functions -----------------------------------------------------

function [xi,ri,xs,rs] = aux_pop(xs,rs,bs)   
    bs = min(bs,size(xs,2));

    % Pop first bs elements from xs.
    idx = 1:bs;
    % Pop random elements.
    % idx = idx(randperm(length(idx)));

    xi = xs(:,idx);
    xs(:,idx) = [];

    ri = rs(:,idx);
    rs(:,idx) = [];
end

function [xis,ris] = aux_split(xi,ri,sens,nSplits,nDims)
    [n,bs] = size(xi);
    % Cannot split more than every dimension.
    nDims = min(n,nDims);
    % Split each input in the batch into nSplits parts; use radius*sens 
    % as the splitting heuristic.
    % 1. Find the input dimension with the largest heuritic.
    [~,sortDims] = sort(abs(sens.*ri),1,'descend');
    dimIds = sortDims(1:nDims,:); 

    splitsIdx = repmat(1:nSplits,1,bs);
    bsIdx = repelem((1:bs)',nSplits);

    dim = dimIds(1,:);
    linIdx = sub2ind([n bs nSplits], ...
        repelem(dim,nSplits),bsIdx(:)',splitsIdx(:)');
    % 2. Split the selected dimension.
    xi_ = xi;
    ri_ = ri;
    % Shift to the lower bound.
    dimIdx = sub2ind([n bs],dim,1:bs);
    xi_(dimIdx) = xi_(dimIdx) - ri(dimIdx);
    % Reduce raidus.
    ri_(dimIdx) = ri_(dimIdx)/nSplits;
   
    xis = repmat(xi_,1,1,nSplits);
    ris = repmat(ri_,1,1,nSplits);
    % Offset the center.
    xis(linIdx(:)) = xis(linIdx(:)) + (2*splitsIdx(:) - 1).*ris(linIdx(:));
    
    % xis = permute(xis,[1 2 4 3]);
    % ris = permute(ris,[1 2 4 3]);
    % Flatten.
    xis = xis(:,:);
    ris = ris(:,:);
    % % Update number of batch entries.
    % bs = size(xis,2);
end

% ------------------------------ END OF CODE ------------------------------
