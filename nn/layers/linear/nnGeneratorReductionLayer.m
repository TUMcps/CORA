classdef nnGeneratorReductionLayer < nnIdentityLayer
% nnGeneratorReductionLayer - reduces the number of generators of a zonotope
%
% Syntax:
%    obj = nnGeneratorReductionLayer(maxGens)
%    obj = nnGeneratorReductionLayer(maxGens, name)
%
% Inputs:
%    maxGens - maximum order
%    name - name of the layer, defaults to type
%
% Outputs:
%    obj - generated object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork

% Authors:       Lukas Koller
% Written:       30-April-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties
    maxGens                % maximum number of generators
end

methods
    % constructor
    function obj = nnGeneratorReductionLayer(maxGens, varargin)
        % parse input
        inputArgsCheck({ ...
            {maxGens, 'att', 'numeric', 'gpuArray'}; ...
        })

        obj.maxGens = maxGens;
    end
end

% evaluate ----------------------------------------------------------------

methods  (Access = {?nnLayer, ?neuralNetwork})
    
    % zonotope batch (for training)
    function [c, G] = evaluateZonotopeBatch(obj, c, G, options)
        % Reduce the generator matrix.
        [G_,I,keepGensIdx,reduceGensIdx] = aux_reduceGirad(obj,G,obj.maxGens);
        [n,q,batchSize] = size(G_);

        % Compute indices for approximation errors in the generator
        % matrix.
        GdIdx = reshape(sub2ind([n q+n batchSize], ...
            repmat(1:n,1,batchSize), ...
            repmat(1:n,1,batchSize), ...
            repelem(1:batchSize,1,n)),[n batchSize]);

        % Append generators for the approximation errors.    
        G = cat(2,G_,zeros(n,n,batchSize,'like',G_));
        % Add approximation errors to the generators.
        G(GdIdx) = I;

        % Store the indices of the reduced generators.
        if options.nn.train.backprop
            obj.backprop.store.keepGensIdx = keepGensIdx;
            obj.backprop.store.reduceGensIdx = reduceGensIdx;
            obj.backprop.store.GdIdx = GdIdx;
        end

        % % Reduce the generator matrix.
        % [G,U] = aux_reducePCA(obj,G);
        % 
        % % Store the indices of the reduced generators.
        % if options.nn.train.backprop
        %     obj.backprop.store.U = U;
        % end
    end

    % backprop ------------------------------------------------------------

    function [gc, gG] = backpropZonotopeBatch(obj, c, G, gc, gG, options)
        % Extract indices of reduced generators.
        keepGensIdx = obj.backprop.store.keepGensIdx;
        reduceGensIdx = obj.backprop.store.reduceGensIdx;
        GdIdx = obj.backprop.store.GdIdx;

        [n,~,batchSize] = size(G);
        Gred = reshape(G(:,reduceGensIdx),n,[],batchSize);

        gG_ = zeros(size(G),'like',G);
        gG_(:,keepGensIdx) = reshape(gG(:,1:size(keepGensIdx,1),:),n,[]);
        gG_(:,reduceGensIdx) = reshape(sign(Gred).*permute(gG(GdIdx),[1 3 2]),n,[]);
        gG = gG_;

        % % Extract indices of reduced generators.
        % U = obj.backprop.store.U;
        % 
        % [n,~,batchSize] = size(G);
        % 
        % UgG = pagemtimes(U,'transpose',gG,'none');
        % % Compute indices for diagonal entries in the generator matrix.
        % dIdx = reshape(sub2ind([n n batchSize], ...
        %     repmat(1:n,1,batchSize), ...
        %     repmat(1:n,1,batchSize), ...
        %     repelem(1:batchSize,1,n)),[n batchSize]);
        % 
        % gG = pagemtimes(U,sign(G).*permute(UgG(dIdx),[1 3 2]));
    end
end

% Auxiliary functions -----------------------------------------------------

methods

    function [G_,I,keepGensIdx,reduceGensIdx] = aux_reduceGirad(obj,G,maxGens)
        % Obtain number of dimensions and generators.
        [n,q,batchSize] = size(G);
        % Compute the length of each generator.
        genLens = reshape(sum(abs(G)),[q batchSize]);
        % Sort the generators by their length.
        [~,idx] = sort(genLens);
        idx = reshape(sub2ind([q batchSize], ...
            idx,repmat(1:batchSize,q,1)),size(idx));
        % Identify the generators to keep.
        keepGensIdx = idx(1:maxGens - n,:);
        G_ = reshape(G(:,keepGensIdx),n,[],batchSize);
        % Identify generators to reduce.
        reduceGensIdx = idx(maxGens - n + 1:end,:);
        Gred = reshape(G(:,reduceGensIdx),n,[],batchSize);
        I = reshape(sum(abs(Gred),2),[n batchSize]);
    end

    function [G_,U] = aux_reducePCA(obj,G)
        % Obtain number of dimensions and generators.
        [n,~,batchSize] = size(G);

        X = pagetranspose(cat(2,G,-G));
        Co = pagemtimes(X,'transpose',X,'none');
        U = zeros([n n batchSize],'like',G);
        for i=1:batchSize
            [U(:,:,i),~,~] = svd(Co(:,:,i));
        end
        % [U,~,~] = pagesvd(Co);
        r = sum(abs(pagemtimes(U,'transpose',G,'none')),2);
        idmat = eye(size(r,1),'like',G);
        G_ = pagemtimes(U,r.*idmat);
    end
end

end

% ------------------------------ END OF CODE ------------------------------
