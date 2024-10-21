function numGen = prepareForZonoBatchEval(nn,x,varargin)
% prepareForZonoBatchEval - perpare neural network for zonotope batch
%   evaluation: In each layer store the active generator ids and identity 
%   matrices to quickly add approximation errors.
%
% Syntax:
%    nn.prepareForZonoBatchEval(x, numInitGens, useApproxError, idxLayer)
%
% Inputs:
%    x - example input; used for size and type
%    options
%    idxLayer - indices of layers that should be evaluated
%
% Outputs:
%    numGen - number of generators used during training
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork/evaluateZonotopeBatch

% Authors:       Lukas Koller
% Written:       07-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

v0 = size(x,1);

% validate parameters
[options,idxLayer] = setDefaultValues({struct,1:length(nn.layers)}, varargin);
options = nnHelper.validateNNoptions(options,true);

% compute total number of generators
if options.nn.train.num_init_gens < inf
    numGen = options.nn.train.num_init_gens;
else
    numGen = v0;
end
% numGen = min(options.nn.train.num_init_gens,v0);
numApproxErr = options.nn.train.num_approx_err;
numNeurons = v0; % store number of neurons of the current layer

for i=idxLayer
    % Extract current layer.
    layeri = nn.layers{i};
    % Handle the sublayers of composite layers separately.
    if isa(layeri,'nnCompositeLayer')
        for j=1:length(layeri.layers)
            % Extract the j-th computation path.
            layerij = layeri.layers{j};
            % Iterate over the layers of the current computation path.
            for k=1:length(layerij)
                % Prepare the current layer and update the number of 
                % generators and neurons.
                [numGen,numNeurons] = aux_prepareLayer(layerij{k}, ...
                    numApproxErr,numGen,numNeurons,options);
            end
        end
    else
        % Prepare the current layer and update the number of generators and
        % neurons.
        [numGen,numNeurons] = aux_prepareLayer(layeri,numApproxErr, ...
            numGen,numNeurons,options);
    end
end

end


% Auxiliary functions -----------------------------------------------------

function [numGen,numNeurons] = aux_prepareLayer(layeri,numApproxErr, ...
    numGen,numNeurons,options)
    % Store number generators of the input
    layeri.backprop.store.genIds = 1:numGen;
    if isa(layeri,'nnGeneratorReductionLayer')
        layeri.maxGens = min(numGen,layeri.maxGens);
        numGen = layeri.maxGens;
    elseif isa(layeri,'nnActivationLayer')
        if numApproxErr > 0 && ~options.nn.interval_center
            % % We store an id-matrix in each activation layer to append 
            % % the approximation errors of image enclosures to the 
            % % generator matrix. This eliminates the need to allocate 
            % % new GPU memory during training. 
            % layerIdMat = eye(numNeurons,'like',x);
            % layeri.backprop.store.idMat = layerIdMat; % store id-matrix
            % additionally, we store the indices of the generators
            % corresponding the approximations error of the activation
            % layer. The activation layer simply (i) multiplies the 
            % id-matrix with the vector containing approximation error 
            % and (ii) copies the new generators in the correct spot in
            % the generator matrix.
            layerNumApproxErr = min(numApproxErr,numNeurons);
            layeri.backprop.store.approxErrGenIds = 1 + ...
                numGen:(numGen + layerNumApproxErr);
            % add a generator for each neuron to store the approx. error
            numGen = numGen + min(numApproxErr,numNeurons);
        else
            % no approximation error are stored
            layeri.backprop.store.approxErrGenIds = [];
        end
    elseif isa(layeri,'nnElementwiseAffineLayer')
        % Does not change the number of neurons.
        return
    else
        % Update the number of neurons.
        [~,numNeurons] = layeri.getNumNeurons();
    end
end

% ------------------------------ END OF CODE ------------------------------
