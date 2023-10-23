function res = test_nn_neuralNetwork_getNumNeurons()
% test_nn_neuralNetwork_getNumNeurons - unit test function for 
%     neuralNetwork/getNumNeurons
%
% Syntax:
%    res = test_nn_neuralNetwork_getNumNeurons()
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

resvec = [];

% empty case
nn = neuralNetwork();
numNeurons = nn.getNumNeurons();
resvec(end+1) = length(numNeurons) == 1 && isnan(numNeurons);

% single layer
nn = neuralNetwork({nnLinearLayer([2 1; 0 -4])});
numNeurons = nn.getNumNeurons();
resvec(end+1) = isequal([2,2],numNeurons);
nn = neuralNetwork({nnSigmoidLayer()});
numNeurons = nn.getNumNeurons();
resvec(end+1) = length(numNeurons) == 2 && all(isnan(numNeurons));

% larger network
nn = neuralNetwork({ ...
    nnLinearLayer([2 2; 0 -4; -1 2]); ...
    nnLeakyReLULayer();
    nnLinearLayer([2 2 0; -4 -1 2]); ...
    nnSigmoidLayer()
});
numNeurons = nn.getNumNeurons();
resvec(end+1) = isequal([2,3,3,2,2],numNeurons);

% cnn
nn = neuralNetwork({nnConv2DLayer([2 2; 2 2],2)});
numNeurons = nn.getNumNeurons();
resvec(end+1) = length(numNeurons) == 2 && all(isnan(numNeurons));
% set input size
nn.setInputSize([10 10 1]);
numNeurons = nn.getNumNeurons();
resvec(end+1) = isequal([100,81],numNeurons);

% gather results
res = all(resvec);


% ------------------------------ END OF CODE ------------------------------
