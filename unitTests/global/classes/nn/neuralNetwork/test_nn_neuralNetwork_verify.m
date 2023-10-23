function [res] = test_nn_neuralNetwork_verify()
% test_nn_neuralNetwork_verify - tests the verify function 
%    
%
% Syntax:
%    res = test_nn_neuralNetwork_verify()
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

resvec = [];

% test examples(TODO: make proper unit test)
resvec(end+1) = example_neuralNetwork_verify_01_unsafe_verified();
resvec(end+1) = example_neuralNetwork_verify_02_unsafe_falsified();
resvec(end+1) = example_neuralNetwork_verify_03_safe_verified();
resvec(end+1) = example_neuralNetwork_verify_04_safe_falsified();

% gather results
res = all(resvec);

end

% ------------------------------ END OF CODE ------------------------------
