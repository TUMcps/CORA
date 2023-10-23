function refinable_layers = getRefinableLayers(obj)
% getRefinableLayers - returns all layers which are refinable for
% adaptive evaluation
%
% Syntax:
%    res = getRefinableLayers(obj)
%
% Inputs:
%    obj - object of class neuralNetwork
%
% Outputs:
%    cell array containing refinable layers
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork/evaluate with 'adaptive'

% Authors:       Tobias Ladner
% Written:       28-March-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

layers = obj.layers;
ix = arrayfun(@(layer_i) layer_i{1}.is_refinable, layers);
refinable_layers = layers(ix);

end

% ------------------------------ END OF CODE ------------------------------
