function res = testnn_neuralNetwork_convertDLToolboxNetwork()
% testnn_neuralNetwork_convertDLToolboxNetwork - tests the conversion 
%    to and from networks from the Matlab DL toolbox
%
% Syntax:
%    res = testnn_neuralNetwork_convertDLToolboxNetwork()
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
% Written:       15-May-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% high tol due to DLT using singles
tol = 1e-6;

% test feed-forward neural network ---

% load network
nn = neuralNetwork.readONNXNetwork('nn-nav-set.onnx');
nn_dlt = nn.convertToDLToolboxNetwork();

% test network
x = ones(nn.neurons_in,1);
y = nn.evaluate(x);
y_dlt = nn_dlt.predict(x')';
assert(all(withinTol(y,y_dlt,tol)))

% test convolutional neural network ---

% load network
nn = neuralNetwork.readONNXNetwork('vnn_verivital_avgpool.onnx',false,'BCSS');
nn_dlt = nn.convertToDLToolboxNetwork();

% test network
x = ones(nn.neurons_in,1);
y = nn.evaluate(x);
y_dlt = nn_dlt.predict(reshape(x,nn.layers{1}.inputSize))';
assert(all(withinTol(y,y_dlt,tol)))

% Test a neural networks from VNN-COMP ------------------------------------

% Reset the random number generator.
rng('default');
% Specify the model paths.
modelpaths = {
    [CORAROOT '/models/Cora/nn/ACASXU_run2a_1_2_batch_2000.onnx'];
    [CORAROOT '/models/Cora/nn/ACASXU_run2a_5_3_batch_2000.onnx'];
};

for i=1:length(modelpaths)
    % Load the neural network as a DLT network.
    nn_dlt = importNetworkFromONNX(modelpaths{i}, ...
        InputDataFormats='BSSC',NameSpace='DLT_CustomLayers');
    % Load the neural network as a CORA network.
    nn = neuralNetwork.readONNXNetwork(modelpaths{i},false,'BSSC');
    % Generate a random input.
    x = rand(nn.neurons_in,1);
    % Compute the output with the CORA network.
    y = nn.evaluate(x);
    % Compute the output with the DLT network.
    y_dlt = nn_dlt.predict(reshape(x,nn.layers{1}.inputSize))';
    % Check if the results are within the tolerance.
    assert(all(withinTol(y,y_dlt,tol)));
end

% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
