function res = example_neuralNetwork_verify_visual()
% example_neuralNetwork_verify_visual - example for the visualaizing the 
%   verification of a neural networks using the function 
%   neuralNetwork/verify.
%
% Syntax:
%    res = example_neuralNetwork_verify_visual()
%
% Inputs:
%    -
%
% Outputs:
%    res - string, verification result 
%       ['VERIFIED','COUNTEREXAMPLE','UNKNOWN']
%
% References:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Lukas Koller
% Written:       20-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Reset the random number generator.
rng('default');

% Generate random neural network.
nn = neuralNetwork.generateRandom(NrInputs=2,NrOutputs=2, ...
    ActivationFun='relu',NrLayers=3,NrHiddenNeurons=4);
nn.layers(end) = [];

% Specify initial set.    
x = [0; 0]; % center
r = [1; 1]; % radius

% Specify unsafe set specification.
% A = [-5 1; 1 1];       
% b = [-1.8; -1.2];  
A = [-1 1];
b = -2.27;
safeSet = false;

% Verbose verification output.
verbose = true;
% Set a timeout of 2s.
timeout = 2;

% Create evaluation options.
options.nn = struct(...
    'use_approx_error',true,...
    'poly_method','bounds',...'bounds','singh'
    'train',struct(...
        'backprop',false,...
        'mini_batch_size',512 ...
    ) ...
);
% Set default training parameters
options = nnHelper.validateNNoptions(options,true);
options.nn.interval_center = false;

% Set the falsification method: {'fgsm','center','zonotack'}.
options.nn.falsification_method = 'zonotack';
% Set the input set refinedment method: {'naive','zonotack'}.
options.nn.refinement_method = 'naive';

% Do verification.
[res,x_,y_] = nn.verify(x,r,A,b,safeSet,options,timeout,verbose,[1:2; 1:2]);

end


% ------------------------------ END OF CODE ------------------------------
