function numParams = getNumParams(obj)
% getNumParams - returns the total number of learnable parameters in a
%   neural network.
%
% Syntax:
%    pattern = getNumParams(obj)
%
% Inputs:
%    obj - object of class neuralNetwork
%
% Outputs:
%    numParams - number of learnable parameters
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Lukas Koller
% Written:       15-December-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Enumerate all layers.
layers = obj.enumerateLayers();

% Initialize the result.
numParams = 0;

% Iterate over all layers.
for i = 1:length(layers)
    % Obtain the i-th layer.
    layeri = layers{i};
    % Obtain the learnable parameters.
    paramNames = layeri.getLearnableParamNames();
    % Iterate over the parameters.
    for j=1:length(paramNames)
        % Add the number of parameters to the total number.
        numParams = numParams + numel(layeri.(paramNames{j}));
    end
end

end

% ------------------------------ END OF CODE ------------------------------
