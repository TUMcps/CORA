function resetGNN(obj)
% resetGNN - resets the graph neural network properties in each layer
%
% Syntax:
%    resetGNN(obj)
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
% Written:       16-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% delete global message passing property
P = nnGCNLayer.message_passing;
P.val = [];

end

% ------------------------------ END OF CODE ------------------------------
