classdef nnActLayerFromHandle < nnActivationLayer
% nnActLayerFromHandle - activation layer for arbitrary activation function
%    Caution: Use wisely. This only applies general reachability analyzes
%    based on what it can infer from the given function handle. Often, you
%    can make exploit certain properties and get tighter results or 
%    speed up computation time.
%
% Syntax:
%    obj = nnActLayerFromHandle(fun, layerid, name)
%
% Inputs:
%    fun - function handle
%    layerid - identifier of layer (used to re-use derivative computations)
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

% Authors:       Tobias Ladner
% Written:       06-October-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties
    fun = @(x) x; % init with something simple (will be updated later)
end

methods
    % constructor
    function obj = nnActLayerFromHandle(fun,varargin)
        narginchk(1,3)
        [layerid,name] = setDefaultValues({[],[]},varargin);
        % call super class constructor
        obj@nnActivationLayer(name)
        % save property
        inputArgsCheck({{fun,'att','function_handle'}})
        obj.fun = fun;
        % recompute derivative (with unique layer name)
        if isempty(layerid)
            layerid = getDefaultName(obj);
        end
        obj.layerid = layerid;
        obj.df = [];
        obj.df = obj.getDf(1);
    end

    function [df_l,df_u] = getDerBounds(obj, l, u)
        df = obj.df(interval(l,u));
        df_l = df.inf; df_u = df.sup;
    end
end

% evaluate ----------------------------------------------------------------

methods (Access = {?nnLayer, ?neuralNetwork})

    % numeric
    function r = evaluateNumeric(obj, input, options)
        r = obj.fun(input);
    end
    
end
end

% ------------------------------ END OF CODE ------------------------------
