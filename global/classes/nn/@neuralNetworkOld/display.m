function display(obj)
% display - displays the properties of the neuralNetworkOld class 
%
% Syntax:  
%    display(obj)
%
% Inputs:
%    obj - neuralNetworkOld object
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Tobias Ladner
% Written:       23-November-2022
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

% convert to new neuralNetwork
nn = neuralNetwork.getFromOldNeuralNetwork(obj);

% disp input if necessary
dispInput(inputname(1))

display(nn);

%------------- END OF CODE --------------
