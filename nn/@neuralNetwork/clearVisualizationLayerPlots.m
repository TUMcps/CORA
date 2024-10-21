function clearVisualizationLayerPlots(obj)
% clearVisualizationLayerPlots - clear figures of visualization layers
%
% Syntax:
%    clearVisualizationLayerPlots(obj)
%
% Inputs:
%    obj - neuralNetwork
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: NeuralNetwork

% Authors:       Lukas Koller, Tobias Ladner
% Written:       24-June-2022
% Last update:   ---
% Last revision: 17-August-2022

% ------------------------------ BEGIN CODE -------------------------------

n = length(obj.layers);

for i = 1:n
    layer_i = obj.layers{i};
    if isa(layer_i, 'NNVisualizationLayer')
        layer_i.clearPlot();
    end
end

end

% ------------------------------ END OF CODE ------------------------------
