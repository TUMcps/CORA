function res = run_cnn_unittest(modelfile)
% run_cnn_unittest - test basic functionality of cnn evaluation.
%    This unit test checks if the cora model has the same output as the 
%    deep learning toolbox and if propagated samples are inside the 
%    propagate set.
%
% Syntax:
%    res = run_cnn_unittest(modelfile)
%
% Inputs:
%    modelfile - string, onnx model file
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
% Written:       02-December-2022
% Last update:   17-January-2023 (reshape)
%                26-February-2023 (added 'PackageName' pair to ONNX import)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% load an image
digitDatasetPath = fullfile(matlabroot,'toolbox','nnet','nndemos', 'nndatasets','DigitDataset');
imds = imageDatastore(digitDatasetPath, 'IncludeSubfolders',true,'LabelSource','foldernames', 'ReadSize', 1);
X = readimage(imds,3333); % read '3'
X = double(X) / 255;

% load models
try
    nn_dlt = importONNXNetwork(modelfile, 'InputDataFormats', 'BCSS', ...
        'OutputDataFormats', 'BC','TargetNetwork', 'dlnetwork', ...
        'PackageName', 'CustomLayers');
catch
    % Some versions of the ONNX toolbox have issues with importing
    % neural networks with custom layers. Don't throw an error due to
    % this bug.
    warning('Unable to import ONNX network. Might be an issue with custom layers.')
    res = true;
    return;
end
nn_cora = neuralNetwork.convertDLToolboxNetwork(nn_dlt.Layers);

% compare with Deep Learning Toolbox ---------------------------------

% make predictions
X_dlt = dlarray(X,'SSC');
pred_dlt = predict(nn_dlt, X_dlt);
pred_dlt = extractdata(pred_dlt);
pred_dlt = double(pred_dlt);

X_cora = reshape(X, [], 1);
pred_cora = nn_cora.evaluate(X_cora);

% compare results: as DLT works with single instead of double,
% we have to be more liberal ...
res = res & compareMatrices(pred_dlt, pred_cora, 1e-5);

% check set-based predicition -----------------------------------------

% construct input set
delta = 0.005; % adjust perturbation
c = X_cora; % use image as center
dim = length(c); 
% reduce the number of noise generators to speed up the unit test
G = delta*eye(dim, dim/8);
E = eye(dim, dim/8);
pZ = polyZonotope(c, G, [], E);

% propagate set
evParams = struct();
evParams.maxpool_type = "project"; 
pred = nn_cora.evaluate(pZ, evParams);

% propagate samples
xs = [pZ.randPoint(50), pZ.randPoint(10, 'extreme')];
ys = nn_cora.evaluate(xs);

% check if all samples are contained
pred = interval(pred);
res = res && all(pred.contains(ys));

end

% ------------------------------ END OF CODE ------------------------------
