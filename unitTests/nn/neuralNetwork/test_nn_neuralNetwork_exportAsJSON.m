function res = test_nn_neuralNetwork_exportAsJSON()
% test_nn_neuralNetwork_exportAsJSON - unit test function for 
%     neuralNetwork/exportAsJSON
%
% Syntax:
%    res = test_nn_neuralNetwork_exportAsJSON()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       10-November-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

try

    % create network
    nn = neuralNetwork({ ...
        nnLinearLayer([2 3; 4 5]); ...
        nnSigmoidLayer(); ...
        nnLinearLayer([-1 5; 2 -3]); ...
        nnSigmoidLayer(); ...
    });
    
    nnJson = nn.exportAsJSON();
    jsondecode(nnJson); % should be valid json
    resvec(end+1) = true;
    
    % empty case
    nn = neuralNetwork();
    nnJson = nn.exportAsJSON();
    jsondecode(nnJson); % should be valid json
    resvec(end+1) = true;
    
    % random case
    nn = neuralNetwork.generateRandom();
    nnJson = nn.exportAsJSON();
    jsondecode(nnJson); % should be valid json
    resvec(end+1) = true;

catch ME
    % should not fail
    resvec(end+1) = false;
end

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
