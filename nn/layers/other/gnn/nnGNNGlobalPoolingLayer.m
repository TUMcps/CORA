classdef nnGNNGlobalPoolingLayer < nnGNNLayer
% nnGNNGlobalPoolingLayer - class for global pooling gnn layer
%
% Syntax:
%    obj = nnGNNGlobalPoolingLayer(type, name)
%
% Inputs:
%    type - type of global pooling, one of {'add', 'mean'}
%    name - name of the layer, defaults to type
%
% Outputs:
%    obj - generated object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tianze Huang, Gerild Pjetri, Tobias Ladner
% Written:       22-December-2022
% Last update:   23-February-2023 (TL, clean up)
% Last revision: --- 

% ------------------------------ BEGIN CODE -------------------------------

properties (Constant)
    is_refinable = false
end

properties
    type
end

methods
    % constructor
    function obj = nnGNNGlobalPoolingLayer(type, varargin)
        % parse input
        narginchk(0,2)
        name = setDefaultValues({[]}, varargin);
        inputArgsCheck({{type, 'str', {'add', 'mean'}}})

        % call super class constructor
        obj@nnGNNLayer(name)
        obj.type = type;
    end

    function outputSize = getOutputSize(obj, inputSize, G)
        if nargin < 3
            nrNodes = 1;
        else
            nrNodes = G.numnodes;
        end  
        outputSize = [inputSize(1)/nrNodes, 1];
    end

    function [nin, nout] = getNumNeurons(obj)
        nin = [];
        nout = [];
    end
end

% evaluate ----------------------------------------------------------------

methods(Access = {?nnLayer, ?neuralNetwork})

    % numeric
    function r = evaluateNumeric(obj, input, options)
        % init
        nrNodes = options.nn.graph.numnodes;

        % construct matrix
        switch obj.type
            case "add"
                M = ones(1,nrNodes);
            case "mean"
                M = ones(1,nrNodes) * 1/nrNodes;
            otherwise
             
        end 

        % propagate
        r = vecleftmtimes(M,input);
    end

    % sensitivity
    function S = evaluateSensitivity(obj, S, x, options)

        % init
        nrNodes = options.nn.graph.numnodes;
        nrFeat = size(x,1)/nrNodes;

        % construct matrix
        switch obj.type
            case "add"
                M = kron(speye(nrFeat),ones(1, nrNodes));
            case "mean"
                M = kron(speye(nrFeat),ones(1, nrNodes));
                M = M * 1/nrNodes;
            otherwise
                throw(CORAerror("CORA:notSupported",sprintf('Unsupported type ''%s''.',obj.type)))
        end 
        S = S * M;
    end

    % zonotope/polyZonotope
    function [c, G, GI, E, id, id_, ind, ind_] = evaluatePolyZonotope(obj, c, G, GI, E, id, id_, ind, ind_, options)
        
        % init
        nrNodes = options.nn.graph.numnodes;

        % construct matrix
        switch obj.type
            case "add"
                M = ones(1,nrNodes);
            case "mean"
                M = ones(1,nrNodes) * 1/nrNodes;
            otherwise
                throw(CORAerror("CORA:notSupported",sprintf('Unsupported type ''%s''.',obj.type)))
        end 

        % propagate
        c = vecleftmtimes(M,c);
        G = vecleftmtimes(M,G);
        GI = vecleftmtimes(M,GI);
    end
end

end

% ------------------------------ END OF CODE ------------------------------
