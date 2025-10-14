function res = testnn_neuralNetwork_residualConnection()
% testnn_neuralNetwork_residualConnection - construct a neural network with
% a residual connection and check its output against the deep learning
% toolbox.
%
% Syntax:
%    res = testnn_neuralNetwork_residualConnection()
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

% Authors:       Lukas Koller
% Written:       27-June-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Specify the number of input, hidden, and output dimensions.
n0 = 10;
nk = 20;
nK = 5;

% Specify input set.
X0 = interval(zeros(n0,1),ones(n0,1));
% Specify number of inputs to check.
N = 1000;

% Define the layers of the deep learning toolbox network.
layers = [
    featureInputLayer(n0,'Name','input')
    fullyConnectedLayer(nk)
    reluLayer
    fullyConnectedLayer(n0,'Name','fc') % Same as input for residual
    additionLayer(2,'Name','add')
    reluLayer
    fullyConnectedLayer(nK)
];
% Create a layer graph.
lgraph = layerGraph(layers);
% Connect skip connection from input directly to add/in2
lgraph = connectLayers(lgraph,'input','add/in2');

% Convert layer graph to a dlnetwork object and initialize the weights
% randomly.
dlt_net = dlnetwork(lgraph);
% Export the network.
filename = 'dlt_res_net.onnx';
exportONNXNetwork(dlt_net,filename);

% Convert to a CORA network.
nn = neuralNetwork.readONNXNetwork(filename,false,'BC','', ...
    'dagnetwork',true);
% Delete the ONNX file.
delete(filename);

% Sample random inputs.
xs = randPoint(X0,N);

% Compute the outputs.
ys_ = dlt_net.predict(dlarray(xs','BC'));
ys = nn.evaluate(xs);

% Check results.
assert(all(withinTol(ys,ys_,1e-6),'all'));

res = true;

end

% Auxiliary functions -----------------------------------------------------

% ------------------------------ END OF CODE ------------------------------
