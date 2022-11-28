classdef (Abstract) nnActivationLayer < nnLayer
% nnActivationLayer - abstract class for non-linear layers
%
% Syntax:
%    obj = nnActivationLayer(name)
%
% Inputs:
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

% Author:       Tobias Ladner
% Written:      28-March-2022
% Last update:  01-April-2022 (moved to class folder)
% Last revision:10-August-2022 (renamed)

%------------- BEGIN CODE --------------

properties
    % function handles
    f function_handle
    df function_handle
    dfs(:, 1) cell = cell(0, 1)
end

methods
    % constructor
    function obj = nnActivationLayer(name)
        % call super class constructor
        obj@nnLayer(name)

        % init function handles
        obj.f = @(x) obj.evaluateNumeric(x);
        obj.df = obj.getDf(1);

        max_df_init = 1;
        obj.dfs = cell(0, 1);
        for i = 1:max_df_init
            obj.dfs{i, 1} = obj.getDf(i);
        end
    end

    function [nin, nout] = getNumNeurons(obj)
        nin = [];
        nout = [];
    end

    % evaluate
    r = evaluateZonotope(obj, Z, evParams)
    [c, G, Grest, expMat, id, id_, ind, ind_] = evaluatePolyZonotope(obj, c, G, Grest, expMat, id, id_, ind, ind_, evParams)
    [c, G, Grest, expMat, id, id_, ind, ind_] = evaluatePolyZonotopeLin(obj, c, G, Grest, expMat, id, id_, ind, ind_, evParams)
    [c, G, Grest, expMat, id, id_, ind, ind_] = evaluatePolyZonotopeQuad(obj, c, G, Grest, expMat, id, id_, ind, ind_, evParams)
    [c, G, Grest, expMat, id, id_, ind, ind_] = evaluatePolyZonotopeCub(obj, c, G, Grest, expMat, id, id_, ind, ind_, evParams)
    r = evaluateTaylm(obj, input, evParams)
    [c, G, C, d, l, u] = evaluateConZonotope(obj, c, G, C, d, l, u, options, evParams)
end

methods (Abstract)
    evaluateZonotopeNeuron(obj, Z)

    evaluatePolyZonotopeNeuronLin(obj, c_i, G_i, Grest_i, expMat, ...
        ind, ind_, bound_approx)
    evaluatePolyZonotopeNeuronQuad(obj, c_i, G_i, Grest_i, expMat, ...
        ind, ind_, bound_approx)
    evaluatePolyZonotopeNeuronCub(obj, c_i, G_i, Grest_i, expMat, ...
        ind, ind_, bound_approx)

    getDf(obj, i)

    getDerBounds(obj, l, u)

    evaluateTaylmNeuron(obj, input, evParams)
    evaluateConZonotopeNeuron(obj, c, G, C, d, l, u, j, options, evParams)
end

methods (Static)
    layer = instantiateFromString(activation)
end

end

%------------- END OF CODE --------------