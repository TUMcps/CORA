function numMPsteps = getNumMessagePassingSteps(obj)
% getNumMessagePassingSteps - returns the number of message passing steps
%
% Syntax:
%    numMPsteps = getNumMessagePassingSteps(obj)
%
% Inputs:
%    obj - object of class (graph) neuralNetwork

%
% Outputs:
%    numMPsteps - number of message passing steps
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       21-March-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

numMPsteps = sum(cellfun(@(layer) isa(layer,'nnGCNLayer') , obj.layers));

end

% ------------------------------ END OF CODE ------------------------------
