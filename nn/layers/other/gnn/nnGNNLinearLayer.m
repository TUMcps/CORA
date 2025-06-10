classdef nnGNNLinearLayer < nnGNNLayer
% nnGNNLinearLayer - class for linear layers
%
% Syntax:
%    obj = nnGNNLinearLayer(W, b, name)
%
% Inputs:
%    W - weight matrix
%    b - bias column vector
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
%                16-April-2024 (TL, compute as matrix variant)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (Constant)
    is_refinable = false
end

properties
    W, b
end

methods
    % constructor
    function obj = nnGNNLinearLayer(W, varargin)
        % parse input
        [b, name] = setDefaultValues({0, []}, varargin);
        inputArgsCheck({ ...
            {W, 'att', 'numeric'}; ...
            {b, 'att', 'numeric'}; ...
        })

        % check dimensions
        if length(b) == 1
            b = b * ones(size(W, 1), 1);
        end
        if ~all(size(b, 1) == size(W, 1))
           throw(CORAerror('CORA:wrongInputInConstructor', ...
               'The dimensions of W and b should match.'));
        end
        if size(b, 2) ~= 1
           throw(CORAerror('CORA:wrongInputInConstructor', ...
               "Second input 'b' should be a column vector."));
        end

        % call super class constructor
        obj@nnGNNLayer(name)

        obj.W = double(W);
        obj.b = double(b);
    end

    function outputSize = getOutputSize(obj, inputSize, graph)
        if nargin < 3
            nrNodes = 1;
        else
            nrNodes = graph.numnodes;
        end  
        nrOutFeatures = size(obj.W,1);
        outputSize = [nrNodes*nrOutFeatures, 1];
    end

    function [nin, nout] = getNumNeurons(obj)
        nin = size(obj.W, 2);
        nout = size(obj.W, 1);
    end
end

% evaluate ----------------------------------------------------------------

methods  (Access = {?nnLayer, ?neuralNetwork})
    
    % numeric
    function r = evaluateNumeric(obj, input, options)
        r = aux_affineMap(obj,input,true);
    end

    % sensitivity
    function S = evaluateSensitivity(obj, S, x, options)
        Wv = kron(obj.W,speye(options.nn.graph.numnodes));
        S = S * Wv;
        S = full(S);
    end

    % zonotope/polyZonotope
    function [c, G, GI, E, id, id_, ind, ind_] = evaluatePolyZonotope(obj, c, G, GI, E, id, id_, ind, ind_, options)
        c = aux_affineMap(obj,c,true);
        G = aux_affineMap(obj,G,false);
        GI = aux_affineMap(obj,GI,false);
    end
end

% Auxiliary functions -----------------------------------------------------

methods(Access=protected)
    function output = aux_affineMap(obj, input, addBias)
            
        % init
        W = obj.W';
        [k,m] = size(W);
        [nk,h] = size(input);
        n = nk/k;

        % reshape
        input = reshape(input,n,k,h);

        % compute linear map on matrix/matrix set
        if isa(input,'interval')
            output = input * W;
        else
            output = pagemtimes(input,W);
        end


        % add bias
        if addBias
            output = output + obj.b';
        end

        % reshape back to vector shape
        output = reshape(output,n*m,h);
    end
end

end

% ------------------------------ END OF CODE ------------------------------
