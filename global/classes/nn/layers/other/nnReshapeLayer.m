classdef nnReshapeLayer < nnLayer
% nnReshapeLayer - class to reshape the input
%    Usually required between convolutional layers and linear layers
%    or to rearrange between column-major (MATLAB standard) vs. row-major
%    (C standard): https://stackoverflow.com/questions/59793724/reshape-and-indexing-in-matlab-and-python
%
% Syntax:
%    obj = nnReshapeLayer(idx_out)
%
% Inputs:
%    idx_out - indices of reshaped output in the right shape
%
% Outputs:
%    obj - generated object
%
% Example:
%   idx_in = reshape(1:100, 10, 10)
%   idx_out = permute(idx_in, [2, 1])
%   layer = nnReshapeLayer(idx_out)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork

% Authors:       Tobias Ladner
% Written:       17-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (Constant)
    is_refinable = false
end

properties
    idx_out
end

methods
    % constructor
    function obj = nnReshapeLayer(idx_out, varargin)
        obj@nnLayer(varargin{:})
        obj.idx_out = idx_out;
    end

    function outputSize = getOutputSize(obj, inputSize)
        outputSize = size(obj.idx_out);
    end

    function [nin, nout] = getNumNeurons(obj)
        if isempty(obj.inputSize)
            nin = [];
            nout = [];
        else
            % we can only compute the number of neurons if the input
            % size was set.
            nin = prod(obj.inputSize);
            outputSize = getOutputSize(obj, obj.inputSize);
            nout = prod(outputSize);
        end
    end
end

% evaluate ----------------------------------------------------------------

methods (Access = {?nnLayer, ?neuralNetwork})
    
    % numeric
    function r = evaluateNumeric(obj, input, evParams)
        r = obj.aux_reshape(input);
    end

    % sensitivity
    function S = evaluateSensitivity(obj, S, x, evParams)
        S = obj.aux_reshape(S')';
    end

    % zonotope/polyZonotope
    function [c, G, GI, E, id, id_, ind, ind_] = evaluatePolyZonotope(obj, c, G, GI, E, id, id_, ind, ind_, evParams)
        c = obj.aux_reshape(c);
        G = obj.aux_reshape(G);
        GI = obj.aux_reshape(GI);
    end

    % taylm
    function r = evaluateTaylm(obj, input, evParams)
        M = eye(prod(obj.inputSize));
        M = obj.aux_reshape(M);
        r = M * input;
    end

    % conZonotope
    function [c, G, C, d, l, u] = evaluateConZonotope(obj, c, G, C, d, l, u, options, evParams)
        c = obj.aux_reshape(c);
        G = obj.aux_reshape(G);
    end
end

% Auxiliary functions -----------------------------------------------------

methods(Access=private)
    function r = aux_reshape(obj, input)
        idx_vec = reshape(obj.idx_out, [], 1);
        r = input(idx_vec, :);
    end
end

end

% ------------------------------ END OF CODE ------------------------------
