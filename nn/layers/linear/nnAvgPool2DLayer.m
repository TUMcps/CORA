classdef nnAvgPool2DLayer < nnConv2DLayer
% nnAvgPool2DLayer - class for average pooling 2D layers, with
% quadratic pooling region
%
% Syntax:
%    obj = nnAvgPool2DLayer(poolSize, padding, stride, dilation, name)
%
% Inputs:
%    poolSize - dimensions of the pooling region
%    padding - padding [left top right bottom]
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
%    [2] Practical Course SoSe '22 - Report Martina Hinz
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nnConv2DLayer

% Authors:       Martina Hinz, Tobias Ladner
% Written:       17-June-2022
% Last update:   01-December-2022 (combine with nnConv2DLayer)
% Last revision: 17-August-2022

% ------------------------------ BEGIN CODE -------------------------------

properties
    poolSize
end

methods
    %constructor
    function obj = nnAvgPool2DLayer(poolSize, varargin)
        narginchk(1,Inf);
        inputArgsCheck({{poolSize, 'att', 'numeric'}})

        % contruct nnConv2DLayer
        W = ones(poolSize) / prod(poolSize); % dummy filter
        b = 0;
        obj@nnConv2DLayer(W, b, varargin{:})

        obj.poolSize = poolSize;

        % overwrite default values for 
        if nargin < 3
            obj.stride = poolSize;
        end
        if nargin < 5
            name = [];
        end
    end

    % compute size of ouput feature map
    function outputSize = getOutputSize(obj, imgSize)
        in_c = imgSize(3);
        p_h = obj.poolSize(1);
        p_w = obj.poolSize(2);
        % Update filter weights.
        obj.W = 1/(p_h*p_w)*reshape(repelem( ...
            reshape(eye(in_c),1,[]),p_h*p_w,1),[p_h p_w in_c in_c]);

        outputSize = getOutputSize@nnConv2DLayer(obj, imgSize);
        outputSize(end) = imgSize(end); % number of channels remain
    end
end

% internal functions ------------------------------------------------------

methods(Access=protected)
    % function bias = aux_getPaddedBias(obj)
    %     outputSize = getOutputSize(obj, obj.inputSize);
    %     bias = zeros(prod(outputSize), 1);
    % end
end

end

% ------------------------------ END OF CODE ------------------------------
