classdef nnCompositeLayer < nnLayer
% nnCompositeLayer - realize multiple parallel computation paths, e.g. res
%    connections
%
% Syntax:
%    obj = nnCompositeLayer(layers,aggregation)
%
% Inputs:
%    layers - cell array with the different parallel computation paths
%    aggregation - 'add' or 'concant'
%
% Outputs:
%    obj - generated object
%
% Example:
%   % A res connection.
%   obj = nnCompositeLayer({{nnLinearLayer(1,0), nnReLULayer}; {}}, 'add')
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork

% Authors:       Lukas Koller
% Written:       27-June-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (Constant)
    is_refinable = false
end

properties
    layers
    aggregation
end

methods
    % constructor
    function obj = nnCompositeLayer(layers, aggregation, varargin)
        obj@nnLayer(varargin{:})
        obj.layers = layers;
        obj.aggregation = aggregation;
    end

    function outputSize = getOutputSize(obj, inputSize)
        % Compute the output size of the first computation path.
        outputSizes = cell(length(obj.layers),1);
        for i=1:length(obj.layers)
            layersi = obj.layers{i};
            outputSizes{i} = inputSize;
            for j=1:length(layersi)
                outputSizes{i} = layersi{j}.computeSizes(outputSizes{i});
            end
        end

        if strcmp(obj.aggregation,'add')
            % Same output size as the individual outputs.
            outputSize = outputSizes{1};
        elseif strcmp(obj.aggregation,'concat')
            % TODO
        else
            throw(CORAerror('CORA:wrongFieldValue', ...
                'nnCompositeLayer.aggregation', ...
                "Only supported values are 'add' and 'concat'!"));
        end
    end

    function [nin, nout] = getNumNeurons(obj)
        if isempty(obj.inputSize)
            nin = [];
            nout = [];
        else
            % We can only compute the number of neurons if the input
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
        r = 0;
        for i=1:length(obj.layers)
            % Compute result for the i-th computation path.
            layersi = obj.layers{i};
            ri = input;
            for j=1:length(layersi)
                ri = layersi{j}.evaluateNumeric(ri, options);
            end
            % Aggreate results.
            if strcmp(obj.aggregation,'add')
                % Add results.
                r = r + ri;
            elseif strcmp(obj.aggregation,'concat')
                % Concatenate results.
                % TODO
            else
                throw(CORAerror('CORA:wrongFieldValue', ...
                    'nnCompositeLayer.aggregation', ...
                    "Only supported values are 'add' and 'concat'!"));
            end
        end
    end

    % sensitivity
    function S = evaluateSensitivity(obj, S, x, options)
        % Retain input sensitivity.
        Sin = S;
        % Initialize output sensitvity.
        S = 0;
        for i=1:length(obj.layers)
            % Compute result for the i-th computation path.
            layersi = obj.layers{i};

            % TODO: remove recomputing of intermediate reults; this is 
            % really inefficient
            ris = cell(length(layersi), 1);
            ris{1} = x;
            for j=1:length(layersi)
                ris{j+1} = layersi{j}.evaluateNumeric(ris{j},options);
            end

            Si = Sin;
            for j=length(layersi):-1:1
                Si = layersi{j}.evaluateSensitivity(Si,ris{j},options);
            end
            % Aggreate results.
            if strcmp(obj.aggregation,'add')
                % Add results.
                S = S + Si;
            elseif strcmp(obj.aggregation,'concat')
                % Concatenate results.
                % TODO
            else
                throw(CORAerror('CORA:wrongFieldValue', ...
                    'nnCompositeLayer.aggregation', ...
                    "Only supported values are 'add' and 'concat'!"));
            end
        end
    end

    % zonotope batch
    function [rc, rG] = evaluateZonotopeBatch(obj, c, G, options)
        rc = 0;
        rG = 0;
        for i=1:length(obj.layers)
            % Initialize output of the i-th computation path. 
            rci = c;
            rGi = G;
            % Compute result for the i-th computation path.
            layersi = obj.layers{i};
            for j=1:length(layersi)
                [rci,rGi] = layersi{j}.evaluateZonotopeBatch(rci,rGi,options);
            end

            % Aggreate results.
            if strcmp(obj.aggregation,'add')
                % Add final results.
                rc = rc + rci;
                if size(rG,2) < size(rGi,2)
                    rGi(:,1:size(rG,2),:) = rGi(:,1:size(rG,2),:) + rG;
                    rG = rGi;
                else
                    rG(:,1:size(rGi,2),:) = rG(:,1:size(rGi,2),:) + rGi;
                end
            elseif strcmp(obj.aggregation,'concat')
                % Concatenate results.
                % TODO
            else
                throw(CORAerror('CORA:wrongFieldValue', ...
                    'nnCompositeLayer.aggregation', ...
                    "Only supported values are 'add' and 'concat'!"));
            end
        end
    end
end

% Auxiliary functions -----------------------------------------------------

end

% ------------------------------ END OF CODE ------------------------------
