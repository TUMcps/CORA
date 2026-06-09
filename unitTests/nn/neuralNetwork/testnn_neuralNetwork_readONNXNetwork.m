function res = testnn_neuralNetwork_readONNXNetwork()
% testnn_neuralNetwork_readONNXNetwork - tests the readONNXNetwork function
%
% Syntax:
%    res = testnn_neuralNetwork_readONNXNetwork()
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

% Authors:       Tobias Ladner, Benedikt Kellner
% Written:       28-November-2022
% Last update:   14-March-2026 (BK, added DLT-vs-CORA roundtrip tests)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% test reading basic network
nn = neuralNetwork.readONNXNetwork('attitude_control_3_64_torch.onnx');

% test verbose output + input/output formats
nn = neuralNetwork.readONNXNetwork('controller_airplane.onnx', true, 'BC', 'BC');

% Reading network with custom layer
nn = neuralNetwork.readONNXNetwork('vnn_verivital_avgpool.onnx', false, 'BCSS');

% high tol due to DLT using single precision internally
tol = 1e-4;

% test ONNX import roundtrip: safenlp (triggers ScalingLayer conversion)
aux_readAndCompare('perturbations_0.onnx', tol);

% test ONNX import roundtrip: MNIST (triggers ScalingLayer conversion)
aux_readAndCompare('mnist-set.onnx', tol);

% gather results
res = true;

end


% Auxiliary functions -----------------------------------------------------

function aux_readAndCompare(onnxFile, tol)
    % load ONNX via CORA (goes through convertDLToolboxNetwork)
    nn = neuralNetwork.readONNXNetwork(onnxFile, false, 'BC');

    % load same ONNX via DLT directly (ground truth)
    dlt_net = importNetworkFromONNX(which(onnxFile), ...
        'InputDataFormats', 'BC', 'OutputDataFormats', 'BC', ...
        'NameSpace', 'DLT_CustomLayers');

    % skip CustomOutputLayer (CORA ignores it during conversion)
    if isa(dlt_net.Layers(end), 'nnet.onnx.layer.CustomOutputLayer')
        args = {'Outputs', dlt_net.Layers(end-1).Name};
    else
        args = {};
    end

    % evaluate and compare
    x = ones(nn.neurons_in, 1);
    y_cora = nn.evaluate(x);
    y_dlt = reshape(double(extractdata( ...
        predict(dlt_net, dlarray(x', 'BC'), args{:}))), [], 1);
    assert(all(withinTol(y_cora, y_dlt, tol), 'all'))
end

% ------------------------------ END OF CODE ------------------------------
