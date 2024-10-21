function resetBounds(obj)
% resetBounds - resets the boudns in each layer
%
% Syntax:
%    resetBounds(obj)
%
% Inputs:
%    obj - object of class neuralNetwork
%
% Outputs:
%    None
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork/reset

% Authors:       Tobias Ladner
% Written:       30-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% reset order
refinable_layers = obj.getRefinableLayers();
for i = 1:length(refinable_layers)
    refinable_layers{i}.l = [];
    refinable_layers{i}.u = [];
end

end

% ------------------------------ END OF CODE ------------------------------
