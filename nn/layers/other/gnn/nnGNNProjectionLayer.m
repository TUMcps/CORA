classdef nnGNNProjectionLayer < nnGNNLayer
% nnGNNProjectionLayer - projection layer to remove irrelevant nodes
%    for node level output
%
% Syntax:
%    obj = nnGNNProjectionLayer(idx_keep,numNodes)
%
% Inputs:
%    nodes_keep - indices of kept nodes
%    numNodes - number of input nodes
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

% Authors:       Tianze Huang, Gerild Pjetri, Tobias Ladner
% Written:       22-December-2022
% Last update:   23-February-2023 (TL, clean up)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (Constant)
    is_refinable = false
end

properties
    P
    nodes_keep
end

methods
    % constructor
    function obj = nnGNNProjectionLayer(nodes_keep,numNodes,varargin)
        % parse input
        [name] = setDefaultValues({[]}, varargin);
        inputArgsCheck({ ...
            {nodes_keep, 'att', {'numeric','logical'},{'vector'}}; ...
            {numNodes, 'att', 'numeric',{'scalar','integer'}}; ...
        })

        % init project matrix
        P = eye(numNodes);
        P = P(nodes_keep,:);

        % call super class constructor
        obj@nnGNNLayer(name)
        obj.P = P;
        obj.nodes_keep = nodes_keep;
    end

    function outputSize = getOutputSize(obj, inputSize)
        featnodes = inputSize(1);
        [nodesIn,nodesOut] = size(obj.P);
        outputSize = [featnodes*nodesOut/nodesIn, 1];
    end

    function [nin, nout] = getNumNeurons(obj)
        nin = [];
        nout = [];
    end
end

% evaluate ----------------------------------------------------------------

methods  (Access = {?nnLayer, ?neuralNetwork})
    
    % numeric
    function r = evaluateNumeric(obj, input, options)
        inputDim = size(input, 1);
        Pfull = aux_computeFullProjectionMatrix(obj, options.nn.graph, inputDim);
        r = Pfull * input;
        r = full(r);
    end

    % sensitivity
    function S = evaluateSensitivity(obj, S, x, options)
        inputDim = size(x, 1);
        Pfull = aux_computeFullProjectionMatrix(obj, options.nn.graph, inputDim);
        S = S * Pfull;
    end

    % zonotope/polyZonotope
    function [c, G, GI, E, id, id_, ind, ind_] = evaluatePolyZonotope(obj, c, G, GI, E, id, id_, ind, ind_, options)
        inputDim = size(c, 1);
        Pfull = aux_computeFullProjectionMatrix(obj, options.nn.graph, inputDim);
        c = Pfull * c;
        G = Pfull * G;
        GI = Pfull * GI;
        
        % make it full again
        c = full(c);
        G = full(G);
        GI = full(GI);
    end
end

% Auxiliary functions -----------------------------------------------------

methods(Access=protected)
    function Pfull = aux_computeFullProjectionMatrix(obj, graph, inputDim)
        % computes vectorized projection matrix
        nrNodes = graph.numnodes;
        nrFeatures = inputDim/nrNodes;

        Pfull = kron(speye(nrFeatures),obj.P);
    end
end

methods
    function updateMessagePassing(obj)
        gcnLayer = nnGCNLayer();

        MP = gcnLayer.message_passing.val;
        if isnumeric(MP)
            if ~isempty(MP)
                % project by remaining nodes
                MP = MP(obj.nodes_keep,:);
                MP = MP(:,obj.nodes_keep);
                gcnLayer.message_passing.val = MP;
            end
        else
            % polyZonotope, do left & right multiplication with projection
            % P * MP * P'
            P = obj.P;
            [numNodesOut,numNodesIn] = size(P); 
            MP = kron(speye(numNodesIn),P) * MP;
            MP = kron(P,speye(numNodesOut)) * MP;
            gcnLayer.message_passing.val = MP;
        end

    end
end

end

% ------------------------------ END OF CODE ------------------------------
