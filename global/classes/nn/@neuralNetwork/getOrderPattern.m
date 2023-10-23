function pattern = getOrderPattern(obj)
% getOrderPattern - returns the current order pattern
%
% Syntax:
%    pattern = getOrderPattern(obj)
%
% Inputs:
%    obj - object of class neuralNetwork
%
% Outputs:
%    pattern - cell array
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork/refine

% Authors:       Tobias Ladner
% Written:       26-August-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

pattern = {};

refinable_layers = obj.getRefinableLayers();
for i = 1:length(refinable_layers)
    pattern{i} = refinable_layers{i}.order';
end

end

% ------------------------------ END OF CODE ------------------------------
