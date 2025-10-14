classdef nnVisualizationLayer < nnLayer
% nnVisualizationLayer - class to visualize intermediate results in
% a network.
%
% Syntax:
%    obj = nnVisualizationLayer(id, visNeuronIds, name)
%
% Inputs:
%    name - name of the layer, defaults to type
%    id - identification number for layer
%    visNeuronIds - (optional) array of ids of neurons to visualize; default: [1 2]
%
% Outputs:
%    obj - generated object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: NeuralNetwork

% Authors:       Lukas Koller, Tobias Ladner
% Written:       23-June-2022
% Last update:   02-October-2023 (TL, generalized plotting)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (Constant)
    is_refinable = false
end

properties
    id, visNeuronIds,
    % store function handles to plotting functions
    visualizeNeuronsPolyZonotope,
    visualizeNeuronsNumeric
end

methods
    % constructor
    function obj = nnVisualizationLayer(id, visNeuronIds, name)
        if nargin < 3
            name = [];
        end
        % call super class constructor
        obj@nnLayer(name)

        if nargin < 2
            visNeuronIds = [1; 2];
        end

        obj.name = name;
        obj.id = id;
        obj.visNeuronIds = visNeuronIds;
        % set plotting function handles
        obj.visualizeNeuronsNumeric = @obj.visualizeNeuronsNumericDefault;
        obj.visualizeNeuronsPolyZonotope = @obj.visualizeNeuronsPolyZonotopeDefault;
    end

    % Return the id for the plot (>10000).
    function fid = getFigureId(obj)
        fid = 10000 + obj.id;
    end

    % clear the plot.
    function r = clearPlot(obj)
        fid = getFigureId(obj);
        figure(fid);
        hold on;
        clf
    end

    % Default plotting function for numerical samples.
    function han = visualizeNeuronsNumericDefault(obj, input)
        obj.setupFigure();
        han = scatter(input(1, :), input(2, :), '.','DisplayName','numeric');
    end

    % Default plotting function for sets
    function han = visualizeSet(obj, S)
        obj.setupFigure();
        han = plot(S, [1, 2],'DisplayName',class(S));
    end

    function f = setupFigure(obj)
        % get figure id
        fid = getFigureId(obj);

        if ~ishandle(fid)
            % set up figure
            f = figure(fid);
            hold on;
            title(obj.name, 'Interpreter', 'none')
            xlabel(['Neuron ', num2str(obj.visNeuronIds(1))])
            ylabel(['Neuron ', num2str(obj.visNeuronIds(2))])
            legend();
        else
            % return existing figure
            f = figure(fid);
        end
    end

    function [nin, nout] = getNumNeurons(obj)
        % returns number of in- and output neurons
        nin = [];
        nout = [];
    end

    function outputSize = getOutputSize(obj, inputSize)
        % returns output size
        outputSize = inputSize;
    end
end

% evaluate ----------------------------------------------------------------

methods (Access = {?nnLayer, ?neuralNetwork})

    % numeric 
    function input = evaluateNumeric(obj, input, options)
        % visualize neurons
        obj.visualizeNeuronsNumeric(input(obj.visNeuronIds, :));
    end

    % interval
    function input = evaluateInterval(obj, input, options)
        % visualize neurons
        obj.visualizeSet(input(obj.visNeuronIds, :));
    end

    % sensitivity
    function S = evaluateSensitivity(obj, S, x, options)
        % return identity
    end

    % zonotope/polyZonotope
    function [c, G, GI, E, id, id_, ind, ind_] = evaluatePolyZonotope(obj, c, G, GI, E, id, id_, ind, ind_, options)
        % visualize neurons
        pZ = polyZonotope(c, G, GI, E, id);
        pZ = project(pZ, obj.visNeuronIds);
        obj.visualizeSet(pZ);
    end
end

end

% ------------------------------ END OF CODE ------------------------------
