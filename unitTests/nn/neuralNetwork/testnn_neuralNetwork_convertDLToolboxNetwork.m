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

% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
