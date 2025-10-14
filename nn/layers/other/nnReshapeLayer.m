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
    function r = evaluateNumeric(obj, input, options)
        r = obj.aux_reshape(input);
    end

    % sensitivity
    function S = evaluateSensitivity(obj, S, options)
        inSize = obj.inputSize;
        S_ = permute(S,[2 1 3]);
        S_ = obj.aux_embed(S_,inSize);
        S = permute(S_,[2 1 3]);

        if options.nn.store_sensitivity
            % Store the gradient (used for the sensitivity computation).
            obj.sensitivity = S;
        end
    end

    % zonotope batch (for training)
    function [c, G] = evaluateZonotopeBatch(obj, c, G, options)
        c = obj.aux_reshape(c);
        G = obj.aux_reshape(G);
    end

    % zonotope/polyZonotope
    function [c, G, GI, E, id, id_, ind, ind_] = evaluatePolyZonotope(obj, c, G, GI, E, id, id_, ind, ind_, options)
        c = obj.aux_reshape(c);
        G = obj.aux_reshape(G);
        GI = obj.aux_reshape(GI);
    end

    % taylm
    function r = evaluateTaylm(obj, input, options)
        M = eye(prod(obj.inputSize));
        M = obj.aux_reshape(M);
        r = M * input;
    end

    % conZonotope
    function [c, G, C, d, l, u] = evaluateConZonotope(obj, c, G, C, d, l, u, options)
        c = obj.aux_reshape(c);
        G = obj.aux_reshape(G);
    end

    % backprop ------------------------------------------------------------

    % numeric
    function grad_in = backpropNumeric(obj, input, grad_out, options, updateWeights)
        inSize = obj.inputSize;
        grad_in = obj.aux_embed(grad_out,inSize);
    end

    % interval batch
    function [gl, gu] = backpropIntervalBatch(obj, l, u, gl, gu, options, updateWeights)
        inSize = obj.inputSize;
        gl = obj.aux_embed(gl,inSize);
        gu = obj.aux_embed(gu,inSize);
    end
    
    % zonotope batch
    function [gc, gG] = backpropZonotopeBatch(obj, c, G, gc, gG, options, updateWeights)
        inSize = obj.inputSize;
        gc = obj.aux_embed(gc,inSize);
        gG = obj.aux_embed(gG,inSize);
    end
end

% Auxiliary functions -----------------------------------------------------

methods(Access=private)
    function r = aux_reshape(obj, input)
        isMatrix = ndims(input) > 2;
        % Obtain the batch size.
        if isMatrix
            [~,q,bSz] = size(input);
            % Reshape input for easier handling.
            input = input(:,:);
        end

        idx_vec = obj.idx_out(:);
        r = input(idx_vec,:);

        if isMatrix
            % Reshape result to original shape.
            r = reshape(r,[],q,bSz);
        end
    end

    function r = aux_embed(obj, input, inSize)
        isMatrix = ndims(input) > 2;
        % Obtain the batch size.
        if isMatrix
            [~,q,bSz_] = size(input);
            % Reshape input for easier handling.
            bSz = q*bSz_;
            input = input(:,:);
        else
            [~,bSz] = size(input);
        end
        % The inverse of reshape; needed for backpropagation.
        idx_vec = obj.idx_out(:);
        r = zeros([prod(inSize) bSz],'like',input);
        r(idx_vec,:) = input;

        if isMatrix
            % Reshape result to original shape.
            r = reshape(r,[],q,bSz_);
        end
    end
end

methods
    function fieldStruct = getFieldStruct(obj)
        fieldStruct = struct;
        fieldStruct.idx_out = obj.idx_out;
    end
end

end

% ------------------------------ END OF CODE ------------------------------
