function [res] = testnn_neuralNetwork_verify()
% testnn_neuralNetwork_verify - tests the verify function 
%    
%
% Syntax:
%    res = testnn_neuralNetwork_verify()
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       01-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Use two acasxu instances from vnn-comp 2023 for testing.
resSafe = example_neuralNetwork_verify_safe(); % Verify the specifiation.
assert(strcmp(resSafe,'VERIFIED'));
resUnsafe = example_neuralNetwork_verify_unsafe(); % Find a counterexample.
assert(strcmp(resUnsafe,'COUNTEREXAMPLE'));

% gather results
res = true;

end

% ------------------------------ END OF CODE ------------------------------
