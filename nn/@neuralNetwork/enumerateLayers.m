function [layersEnum,ancIdx,predIdx,succIdx] = enumerateLayers(nn)
% enumerateLayers - enumerate the layer of a neural network.
%
% Syntax:
%    layersEnum = nn.enumerateLayers()
%
% Inputs:
%    nn - object of class neuralNetwork
%
% Outputs:
%    layersEnum - cell-array of layers
%    ancIdx - index in the top-level ancestor of the layer
%    predIdx - index of predecessor in computation graph
%    succIdx - index of successor in computation graph
%
% References:---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork

% Authors:       Lukas Koller
% Written:       17-April-2025
% Last update:   ---
% Last revision: ---  

% ------------------------------ BEGIN CODE -------------------------------

% Obtain the layers of the neural network.
layers = nn.layers;
% Initialize ancestor indices.
ancIdx = 1:length(layers);

% Enumerate all layers of the neural network.
[layersEnum,ancIdx,predIdx,succIdx] = aux_enumerateLayers(layers,ancIdx);

end


% Auxiliary functions -----------------------------------------------------

function [layersEnum,ancIdx,predIdx,succIdx] = ...
    aux_enumerateLayers(layers,ancIdx)
    % Initialize result.
    layersEnum = {};
    ancIdxIds = []; % Indices into ancIdx.
    predIdx = []; % Indices to predecessors.
    succIdx = []; % Indices to successors.
    % Recursive iteration over the layers.
    for i=1:length(layers)
        % Obtain i-th layer.
        layeri = layers{i};
        % Append the layers of 
        if isa(layeri,'nnCompositeLayer')
            % Store index for predecessor layer.
            predi = length(layersEnum);
            % Enumerate layers of computation paths.
            for j=1:length(layeri.layers)
                % Obtain the top-level layers of the j-th computation path.
                layersij = layeri.layers{j};
                if isempty(layersij)
                    % This is a residual connection. There are now layers
                    % here.
                    continue;
                end
                % All layers in this computation path have the same
                % ancestor.
                ancIdxIdij = repelem(i,1,length(layersij));
                % Enumerate layers of the j-th computation path.
                [layersijEnum,ancIdxIdij,predIdxij,succIdxij] = ...
                    aux_enumerateLayers(layeri.layers{j},ancIdxIdij);
                % Offset the predecessor and successor indices.
                predIdxij = predIdxij + length(layersEnum);
                predIdxij(1) = predi;
                succIdxij = succIdxij + length(layersEnum);
                succIdxij(end) = NaN;
                % Append layers.
                layersEnum = [layersEnum layersijEnum];
                % Replicate ancestor index for all new layers.
                ancIdxIds = [ancIdxIds ancIdxIdij];
                % Append the predecessor and successor indices.
                predIdx = [predIdx predIdxij];
                succIdx = [succIdx succIdxij];
            end
            % Store index of the successor layer.
            succIdx(isnan(succIdx)) = length(layersEnum)+1;
        else
            % Append layer.
            layersEnum{end+1} = layeri;
            ancIdxIds(end+1) = i;
            predIdx(end+1) = length(layersEnum)-1;
            succIdx(end+1) = length(layersEnum)+1;
        end
    end
    % Construct ancestor indices.
    ancIdx = ancIdx(ancIdxIds);
end

% ------------------------------ END OF CODE ------------------------------
