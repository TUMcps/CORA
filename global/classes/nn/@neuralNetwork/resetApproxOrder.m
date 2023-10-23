function resetApproxOrder(obj)
% resetApproxOrder - resets the approximation order in each layer
%
% Syntax:
%    resetApproxOrder(obj)
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
% Written:       25-April-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% reset order
refinable_layers = obj.getRefinableLayers();
for i = 1:length(refinable_layers)
    refinable_layers{i}.order = 1;
end

end

% ------------------------------ END OF CODE ------------------------------
