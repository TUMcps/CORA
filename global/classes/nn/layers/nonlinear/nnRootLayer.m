classdef nnRootLayer < nnActivationLayer
% nnRootLayer - class for square root layers
%
% Syntax:
%    obj = nnRootLayer(name)
%
% Inputs:
%    name - name of the layer, defaults to type
%
% Outputs:
%    obj - generated object
%
% References:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork

% Authors:       Tobias Ladner
% Written:       04-July-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

methods
    % constructor
    function obj = nnRootLayer(name)
        if nargin < 2
            name = [];
        end
        % call super class constructor
        obj@nnActivationLayer(name)
    end

    function df_i = getDf(obj, i)
        if i == 0
            df_i = obj.f;
        else
            df_i1 = obj.getDf(i-1);
            df_i = @(x) 1 / 2 .* df_i1(x);
        end
    end

    function [df_l, df_u] = getDerBounds(obj, l, u)
        df_l = obj.df(l);
        df_u = obj.df(u);
    end
end

% evaluate ----------------------------------------------------------------

methods (Access = {?nnLayer, ?neuralNetwork})
    % numeric
    function r = evaluateNumeric(obj, input, evParams)
        r = sqrt(input);
    end
end
end

% ------------------------------ END OF CODE ------------------------------
