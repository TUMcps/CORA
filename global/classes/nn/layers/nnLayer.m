classdef (Abstract) nnLayer < matlab.mixin.Copyable
% nnLayer - abstract class for nn layers
%
% Syntax:
%    nnLayer
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

% Authors:       Tobias Ladner
% Written:       28-March-2022
% Last update:   15-February-2023 (re-organized, simplified)
%                24-February-2023 (added evParams to all evaluate func)
% Last revision: 10-August-2022 (renamed)

% ------------------------------ BEGIN CODE -------------------------------

properties
    name = [];          % name of the layer

    % stores the size of the input. Needed to propagate images.
    inputSize = [];

    sensitivity = [];   % sensitivity of input
end

methods
    function obj = nnLayer(name)
        if nargin < 1 || isempty(name)
            name = obj.getDefaultName();
        end

        obj.name = name;
        obj.inputSize = [];
    end

    function outputSize = computeSizes(obj, inputSize)
        obj.inputSize = inputSize;
        outputSize = getOutputSize(obj, inputSize);
    end

    function name = getDefaultName(obj)
        name = sprintf("%s_%i", class(obj), nnLayer.getCount());
    end

    function str = getLayerInfo(obj)
        str = sprintf('%-30s %-30s', class(obj), obj.name);
        if ~isempty(obj.inputSize)
            input = obj.inputSize;
            input = join(string(input), 'x');
            output = obj.getOutputSize(obj.inputSize);
            output = join(string(output), 'x');

            str = str + sprintf(" %+10s -> %-10s", input, output);
        else
            str = str + sprintf("\t(Input size not set)");
        end

        % additional information
        if isa(obj, 'nnGNNGlobalPoolingLayer')
            str = str + "(pooling across nodes)";
        elseif isa(obj, 'nnGNNLayer')
            str = str + "(per node)";
        end
    end

    function r = evaluate(obj, input, varargin)
        % wrapper to propagate a single layer
        nn = neuralNetwork({obj});
        r = nn.evaluate(input, varargin{:});
    end
end

% evaluate ----------------------------------------------------------------

methods(Abstract, Access = {?nnLayer, ?neuralNetwork})
    % numeric
    r = evaluateNumeric(obj, input, evParams)
end

methods (Access = {?nnLayer, ?neuralNetwork})
    % evaluation functions for each input set
    % will be overwritten in subclasses if supported by that layer

    % sensitivity
    function S = evaluateSensitivity(obj, S, x, evParams)
        throw(CORAerror('CORA:nnLayerNotSupported', obj, 'sensitivity'))
    end

    % interval
    function bounds = evaluateInterval(obj, bounds, evParams)
        bounds = obj.evaluateNumeric(bounds, evParams);
    end

    % zonotope/polyZonotope
    function [c, G, GI, E, id, id_, ind, ind_] = evaluatePolyZonotope(obj, c, G, GI, E, id, id_, ind, ind_, evParams)
        throw(CORAerror('CORA:nnLayerNotSupported', obj, 'polyZonotope'))
    end

    % taylm
    function r = evaluateTaylm(obj, input, evParams)
        throw(CORAerror('CORA:nnLayerNotSupported', obj, 'taylm'))
    end

    % conZonotope
    function [c, G, C, d, l, u] = evaluateConZonotope(obj, c, G, C, d, l, u, options, evParams)
        throw(CORAerror('CORA:nnLayerNotSupported', obj, 'conZonotope'))
    end
end

% Auxiliary functions -----------------------------------------------------

methods (Abstract)
    [nin, nout] = getNumNeurons(obj)
    outputSize = getOutputSize(obj, inputSize)
end

methods (Access=protected)
    function checkInputSize(obj)
        % checks if input size is set
        if isempty(obj.inputSize)
            throw(CORAerror("CORA:notSupported", 'Input size has to be set to compute CNNs. See neuralNetwork/setInputSize.'))
        end
    end
end

methods (Static, Access = protected)
    function res = getCount()
        % used to give each layer a unique name
        persistent count;
        if isempty(count)
            count = 0;
        end

        count = count + 1;
        res = count;
    end
end

methods (Access = protected)
    % overwrite copy function
    function cp = copyElement(obj)
        cp = copyElement@matlab.mixin.Copyable(obj);
        cp.name = sprintf("%s_copy", obj.name);
    end
end

end

% ------------------------------ END OF CODE ------------------------------
