function gnn_red = reduceGNNForNode(obj,n0,G)
% reduceGNNForNode - removes all nodes that are not relevant to predict the
%   node level output for the given node 
%
% Syntax:
%    gnn_red = reduceGNNForNode(obj,G,n0)
%
% Inputs:
%    obj - object of class (graph) neuralNetwork
%    n0 - numeric, node of interest
%    G - graph
%
% Outputs:
%    gnn_red - cell array
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       21-March-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
inputArgsCheck({ ...
   {obj,'att','neuralNetwork'}; ...
   {n0,'att','numeric',{'scalar','integer'}}; ...
   {G,'att','graph'}; ...
})

numNodes = G.numnodes;
if n0 > numNodes
    throw(CORAerror('CORA:wrongValue','third','n0 has to be part of G'))
end

% init
numMPsteps = obj.getNumMessagePassingSteps();
if numMPsteps == 0
    % no gnn, return
    gnn_red = obj;
end

% init new layers
layers = obj.layers;
layers_red = cell(numel(layers)+numMPsteps+1,1);
cnt = 1;

% compute distance of all neighbors
[remNeighbors,distVec] = nearest(G,n0,G.numnodes);
remNeighbors = [1;remNeighbors];
distVec = [0;distVec];

% initial node reduction
remMPsteps = numMPsteps;
idx_keep = remNeighbors(distVec <= remMPsteps+1);
layers_red{cnt} = nnGNNProjectionLayer(idx_keep,G.numnodes);
remNeighbors = 1:numel(idx_keep);
distVec = distVec(1:numel(idx_keep));
cnt = cnt + 1;

% iterate over network
for k=1:numel(layers)
    layer_k = layers{k};

    % transfer to new layers
    layers_red{cnt} = layer_k.copy();
    cnt = cnt+1;

    if isa(layer_k, 'nnGCNLayer')
        remMPsteps = remMPsteps -1;

        % reduce nodes
        if remMPsteps == 0
            % last MP step computed, only keep node of interest
            idx_keep = remNeighbors(distVec == 0);
        else
            % keep all neighbors in MP + 1 due to normalization
            idx_keep = remNeighbors(distVec <= remMPsteps+1);
        end
        layers_red{cnt} = nnGNNProjectionLayer(idx_keep,numel(remNeighbors));
        remNeighbors = 1:numel(1:numel(idx_keep));
        distVec = distVec(1:numel(idx_keep));
        cnt = cnt + 1;
    end

end

layers_red = layers_red(1:(cnt-1));

gnn_red = neuralNetwork(layers_red);

end

% ------------------------------ END OF CODE ------------------------------
