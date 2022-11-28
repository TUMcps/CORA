function r = evaluateTaylm(obj, input, evParams)
% evaluateTaylm - evaluates the activation layer on a taylor model
%
% Syntax:
%    r = evaluateTaylm(obj, input, evParams)
%
% Inputs:
%    input - taylor model
%    evParams - parameter for NN evaluation
%
% Outputs:
%    r - updated taylor model
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nnLayer

% Author:        Niklas Kochdumper, Tobias Ladner
% Written:       17-September-2021
% Last update:   ---
% Last revision: 05-April-2022 (TL)

%------------- BEGIN CODE --------------

% loop over all neurons in the current layer
r = input;
for j = 1:size(input, 1)
    r(j) = obj.evaluateTaylmNeuron(r(j), evParams);
end

%------------- END OF CODE --------------