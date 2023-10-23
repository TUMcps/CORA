classdef nnConv2DLayer < nnLayer
% nnConv2DLayer - class for convolutional 2D layers
%
% Syntax:
%    obj = nnConv2DLayer(W, b, padding, stride, dilation, name)
%
% Inputs:
%    W - weight matrix (4-D single)
%        all the filters are stored in the weight-matrix:
%        (kernel_height, kernel_width, in_channels, num_filters)
%    b - bias column vector
%    padding - zero padding [left top right bottom]
%    stride - step size per dimension
%    dilation - step size per dimension
%    name - name of the layer, defaults to type
%
% Outputs:
%    obj - generated object
%
% References:
%    [1] T. Gehr, et al. "AI2: Safety and Robustness Certification of
%        Neural Networks with Abstract Interpretation," 2018
%    [2] Practical Course SoSe '22 - Report Lukas Koller
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Lukas Koller, Tobias Ladner
% Written:       04-June-2022
% Last update:   01-December-2022 (combine with nnAvgPool2DLayer)
%                17-January-2023 (TL, optimizations)
% Last revision: 17-August-2022
%                02-December-2022 (clean up)

% ------------------------------ BEGIN CODE -------------------------------

properties (Constant)
    is_refinable = false
end

properties
    W, b, stride, padding, dilation
end

methods
    % constructor
    function obj = nnConv2DLayer(varargin)
        if nargin > 6
            throw(CORAerror('CORA:tooManyInputArgs', 6))
        end
        
        % 1. parse input arguments: varargin -> vars
        [W, b, padding, stride, dilation, name] = aux_parseInputArgs(varargin{:});
        
        % 2. check correctness of input arguments
        aux_checkInputArgs(W, b, padding, stride, dilation, name)

        % 3. call super class constructor
        obj@nnLayer(name)

        % 4. assign properties
        obj.W = W;
        obj.b = b;
        obj.padding = padding;
        obj.stride = stride;
        obj.dilation = dilation;
    end

    function [nin, nout] = getNumNeurons(obj)
        if isempty(obj.inputSize)
            nin = [];
            nout = [];
        else
            % we can only compute the number of neurons if the input
            % size was set.
            nin = obj.inputSize(1) * obj.inputSize(2) * obj.inputSize(3);
            outputSize = getOutputSize(obj, obj.inputSize);
            nout = prod(outputSize);
        end
    end

    % compute size of ouput feature map
    function outputSize = getOutputSize(obj, imgSize)
        in_h = imgSize(1);
        in_w = imgSize(2);
        [f_h, f_w] = aux_getFilterSize(obj);
        % padding [left,top,right,bottom]
        pad_l = obj.padding(1);
        pad_t = obj.padding(2);
        pad_r = obj.padding(3);
        pad_b = obj.padding(4);

        out_h = floor((in_h - f_h + pad_t + pad_b)/obj.stride(1)) + 1;
        out_w = floor((in_w - f_w + pad_l + pad_r)/obj.stride(2)) + 1;
        out_c = size(obj.W, 4);
        outputSize = [out_h, out_w, out_c];
    end
end

% evaluate ------------------------------------------------------------

methods  (Access = {?nnLayer, ?neuralNetwork})
    % All dimensions (h,w,c) are flattened into a vector

    % numeric
    function r = evaluateNumeric(obj, input, evParams)
        % compute weight and bias
        Wff = aux_computeWeightMatrix(obj);
        bias = aux_getPaddedBias(obj);

        % simulate using linear layer
        linl = nnLinearLayer(Wff, bias);
        r = linl.evaluateNumeric(input, evParams);
    end

    % sensitivity
    function S = evaluateSensitivity(obj, S, x, evParams)
        % compute weight and biasinput
        Wff = aux_computeWeightMatrix(obj);
        bias = aux_getPaddedBias(obj);

        % simulate using linear layer
        linl = nnLinearLayer(Wff, bias);
        S = linl.evaluateSensitivity(S, x, evParams);
    end

    % zonotope/polyZonotope
    function [c, G, GI, E, id, id_, ind, ind_] = evaluatePolyZonotope(obj, c, G, GI, E, id, id_, ind, ind_, evParams)
        % compute weight and bias
        Wff = aux_computeWeightMatrix(obj);
        bias = aux_getPaddedBias(obj);

        % simulate using linear layer
        linl = nnLinearLayer(Wff, bias);
        [c, G, GI, E, id, id_, ind, ind_] = linl.evaluatePolyZonotope(c, G, GI, E, id, id_, ind, ind_, evParams);
    end

    % taylm
    function r = evaluateTaylm(obj, input, evParams)
        % compute weight and bias
        Wff = aux_computeWeightMatrix(obj);
        bias = aux_getPaddedBias(obj);

        % simulate using linear layer
        linl = nnLinearLayer(Wff, bias);
        r = linl.evaluateTaylm(obj, input, evParams);
    end

    % conZonotope
    function [c, G, C, d, l, u] = evaluateConZonotope(obj, c, G, C, d, l, u, options, evParams)
        % compute weight and bias
        Wff = aux_computeWeightMatrix(obj);
        bias = aux_getPaddedBias(obj);

        % simulate using linear layer
        linl = nnLinearLayer(Wff, bias);
        [c, G, C, d, l, u] = linl.evaluateConZonotope(obj, c, G, C, d, l, u, options, evParams);
    end
end

% Auxiliary functions ------------------------------------------------------

methods
    function layer = convert2nnLinearLayer(obj)
        % as convolutional layers are just fancy linear layers,
        % they can be converted to them easily.

        % compute weight and bias
        Wff = aux_computeWeightMatrix(obj);
        bias = aux_getPaddedBias(obj);

        layer = nnLinearLayer(Wff, bias, sprintf("%s_linear", obj.name));
    end
end

methods (Access = protected)
    function [f_h, f_w] = aux_getFilterSize(obj)
        % return the size of the filter kernels.
        k_h = size(obj.W, 1);
        k_w = size(obj.W, 2);
        d_h = obj.dilation(1);
        d_w = obj.dilation(2);

        f_h = k_h + (k_h - 1) * (d_h - 1);
        f_w = k_w + (k_w - 1) * (d_w - 1);
    end

    % Compute linear filter-matrix for a single 2D-filter to express a
    % convolution as matrix-vector-multiplication.
    function Wf = aux_computeWeightMatrixForFilter(obj, filter, imgSize)
        % see [1], [2]
        % Only compute the filter-matrix for a single 2D filter acting
        % on only one channel.

        % compute output size
        [f_h, f_w] = size(filter);
        outputSize = getOutputSize(obj, obj.inputSize);
        img_h_out = outputSize(1);
        img_w_out = outputSize(2);

        % compute input size incl. padding
        img_h_in = imgSize(1);
        img_w_in = imgSize(2);        
        % padding [left,top,right,bottom]
        pad_l = obj.padding(1);
        pad_t = obj.padding(2);
        pad_r = obj.padding(3);
        pad_b = obj.padding(4);
        % pad image with zeros
        img_h_in = img_h_in + pad_t + pad_b;
        img_w_in = img_w_in + pad_l + pad_r;

        % stride values
        stride_h = obj.stride(1);
        stride_w = obj.stride(2);

        % turn kernel into vector, padded with 0's
        lin_filter_len = f_h * f_w + (f_h - 1) * (img_h_in - f_h);
        lin_filter = [filter; zeros(img_h_in-f_h, f_w)];
        lin_filter = reshape(lin_filter, [], 1)';
        lin_filter = lin_filter(1:lin_filter_len); % remove trailing zeros

        % init linear weight matrix
        Wf = zeros(img_h_out*img_w_out, img_h_in*img_w_in);
        for i = 1:img_w_out
            k = 1 + (i - 1) * stride_w * img_h_in;
            for j = 1:img_h_out
                % for each row of the weight matrix insert the
                % filter-vector at the correct position.
                idx_out = (i - 1) * img_h_out + j;
                kend = k + lin_filter_len - 1;
                idx_in = k:kend;
                Wf(idx_out, idx_in) = lin_filter;
                k = k + stride_h;
            end
        end

        % remove padded area
        isPad = false(img_h_in, img_w_in);
        isPad(1:pad_t, :) = true;
        isPad(end+1-pad_b:end, :) = true;
        isPad(:, 1:pad_l) = true;
        isPad(:, end+1-pad_r:end) = true;
        isPad = reshape(isPad, 1, []);
        Wf(:, isPad) = [];
    end

    % Computed dilated filter.
    function f = aux_computeDilatedFilter(obj, filter)
        % get filter size
        [k_h, k_w] = size(filter);
        % get dilation factors
        d_h = obj.dilation(1);
        d_w = obj.dilation(2);
        % compute dilated filter size
        df_h = k_h + (k_h - 1) * (d_h - 1);
        df_w = k_w + (k_w - 1) * (d_w - 1);
        f = zeros(df_h, df_w);
        f(1:d_h:end, 1:d_w:end) = filter;
    end

    % Compute weight matrix to express convolutions as matrix-vector multiplication.
    function Wff = aux_computeWeightMatrix(obj)
        checkInputSize(obj)

        % Here we compute the entire weight matrix to express the
        % convolution as a matrix vector multiplication.
        img_h_in = obj.inputSize(1);
        img_w_in = obj.inputSize(2);
        if length(obj.inputSize) < 3 
            c_in = 1;
        else
            c_in = obj.inputSize(3);
        end

        % compute output size
        outputSize = getOutputSize(obj, obj.inputSize);
        img_h_out = outputSize(1);
        img_w_out = outputSize(2);
        if length(outputSize) < 3 
            c_out = 1;
        else
            c_out = outputSize(3);
        end

        % init weight matrix
        Wff = zeros(img_h_out*img_w_out*c_out, img_h_in*img_w_in*c_in);
        for k = 1:c_out % number of filters/out-channels
            for i = 1:c_in % number of in-channels
                if isa(obj, 'nnAvgPool2DLayer')
                    if k == i
                        filter = obj.aux_computeDilatedFilter(obj.W);
                    else
                        % avgPool only within the same channel
                        continue;
                    end
                else
                    % convolve the i-th input channel with the k-th channel
                    filter = obj.aux_computeDilatedFilter(obj.W(:, :, i, k));
                end
                % compute weight matrix for ith filter
                Wf = aux_computeWeightMatrixForFilter(obj, filter, obj.inputSize);
                % compute indices for the ith-filter-weight matrix
                idx1 = (k - 1) * img_h_out * img_w_out + 1;
                idx2 = (i - 1) * img_h_in * img_w_in + 1;
                Wff(idx1:(idx1 + img_h_out * img_w_out - 1), idx2:(idx2 + img_h_in * img_w_in - 1)) = Wf;
            end
        end
    end

    % Pad bias s.t. convolution can be simulated by linear layer.
    function bias = aux_getPaddedBias(obj)
        % We turn the bias into a vector.

        % compute output size
        outputSize = getOutputSize(obj, obj.inputSize);
        img_h_out = outputSize(1);
        img_w_out = outputSize(2);

        % expand the bias vector to output size
        bias = repmat(obj.b', img_h_out*img_w_out, 1);
        bias = reshape(bias, [], 1);
    end
end

end


% Auxiliary functions -----------------------------------------------------

function [W, b, padding, stride, dilation, name] = aux_parseInputArgs(varargin)
    % validate input
    [W, b, padding, stride, dilation, name] = setDefaultValues( ...
        {1, 0, [0 0 0 0], [1 1], [1 1], []}, varargin);
end

function aux_checkInputArgs(W, b, padding, stride, dilation, name)
    % check input types
    inputArgsCheck({ ...
        {W, 'att', 'numeric'}, ...
        {b, 'att', 'numeric'}, ...
        {padding, 'att', 'numeric'}, ...
        {stride, 'att', 'numeric'}, ...
        {dilation, 'att', 'numeric'}, ...
        % name is checked in nnLayer
    })

    % check dimensions
    if length(size(W)) > 4
        throw(CORAerror('CORA:wrongInputInConstructor','Weight matrix has wrong dimensions.'))
    end
    if size(W,4) ~= length(b)
        throw(CORAerror("CORA:wrongInputInConstructor",'Weight matrix and bias dimensions don''t match.'))
    end
    if length(padding) ~= 4
        throw(CORAerror('CORA:wrongInputInConstructor', 'Padding must be an array with 4 entries: [left top right bottom]'))
    end
    if length(stride) ~= 2
        throw(CORAerror('CORA:wrongInputInConstructor', 'Stride must be an array with 2 entries: [dim1, dim2]'))
    end
    if length(dilation) ~= 2
        throw(CORAerror('CORA:wrongInputInConstructor', 'Dilation must be an array with 2 entries: [dim1, dim2]'))
    end

    
end

% ------------------------------ END OF CODE ------------------------------
