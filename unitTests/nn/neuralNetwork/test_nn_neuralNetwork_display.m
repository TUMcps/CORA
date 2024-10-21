function res = test_nn_neuralNetwork_display()
% test_nn_neuralNetwork_display - unit test function for 
%     neuralNetwork/display
%
% Syntax:
%    res = test_nn_neuralNetwork_display()
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
% See also: neuralNetwork/evaluate

% Authors:       Tobias Ladner
% Written:       02-October-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty case
nn = neuralNetwork();
display(nn);

% single layer
nn = neuralNetwork({nnLinearLayer([2 1; 0 -4])});
display(nn);
nn = neuralNetwork({nnSigmoidLayer()});
display(nn);

% larger network
nn = neuralNetwork({ ...
    nnLinearLayer([2 2; 0 -4; -1 2]); ...
    nnLeakyReLULayer();
    nnLinearLayer([2 2 0; -4 -1 2]); ...
    nnSigmoidLayer()
});
display(nn);

% cnn
nn = neuralNetwork({nnConv2DLayer(2,2)});
display(nn);
% set input size
nn.setInputSize([10 10 1]);
display(nn)

res = true;


% ------------------------------ END OF CODE ------------------------------
