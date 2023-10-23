function reset(obj)
% resetBounds - resets the neural network by deleting all values in all
%    layers used for internal computations
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
% See also: neuralNetwork

% Authors:       Tobias Ladner
% Written:       30-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

obj.resetApproxOrder();
obj.resetBounds();

end

% ------------------------------ END OF CODE ------------------------------
