function res = example_neuralNetwork_cnn()
% example_neuralNetwork_cnn - example for the verification of a 
%    convolutional neural network from verital benchmark from VNN'21 [1].
%
% Syntax:
%    res = example_neuralNetwork_cnn()
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean
%
% References:
%    [1] Bak, Stanley, et al. "The second international verification of 
%    neural networks competition (vnn-comp 2021): Summary and results." 
%    arXiv preprint arXiv:2109.00498 (2021).
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       02-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% load an image
digitDatasetPath = fullfile(matlabroot,'toolbox','nnet','nndemos', 'nndatasets','DigitDataset');
imds = imageDatastore(digitDatasetPath, 'IncludeSubfolders',true,'LabelSource','foldernames', 'ReadSize', 1);
X = readimage(imds,3334); % read a '3' image
X = double(X) / 255;

% plot image
figure;
subplot(1, 2, 1);
imshow(X,'InitialMagnification','fit')
title("Input Image")

% reshape to column vector
X = reshape(X, [], 1);

% load model
verbose = false;
nn_cora = neuralNetwork.readONNXNetwork('vnn_verivital_avgpool.onnx', verbose, 'BCSS');
display(nn_cora)

% make predictions
pred_cora = nn_cora.evaluate(X);

disp("Sanity Check:")
disp("    (Label | Prediction)")
disp([(0:9)', pred_cora])
[~, label] = max(pred_cora);
fprintf("Most likely label is: %d\n", label-1);

% check set-based prediction -----------------------------------------

disp("Making a set-based prediction...")

% apply noise to the image
noise = 0.005;

% construct input set
c = X; % use image as center
dim = length(c);
G = noise*eye(dim);
E = eye(dim);
pZ = polyZonotope(c, G, [], E);

% propagate set
evParams = struct();
evParams.maxpool_type = "project"; 
pred = nn_cora.evaluate(pZ, evParams);
pred = interval(pred);

% check results
label_pred = project(pred, label);
other_pred = project(pred, setdiff(1:10, label));

% image + noise is verified if lower bound of the true label is larger
% than the upper bound of all other labels.
isVerified = all(other_pred.sup < label_pred.inf);

if isVerified
    fprintf("VERIFIED with noise=%d.\n", noise)
else
    fprintf("Unable to verify image with noise=%d.\n", noise)
end

% propagate samples
xs = [pZ.randPoint(50), pZ.randPoint(10, 'extreme')];
ys = nn_cora.evaluate(xs);

% plot results ------------------------------------------------------------

subplot(1, 2, 2); hold on;

% plot set prediction
bounds = interval(pred); %,'split');
sups = supremum(bounds);
infs = infimum(bounds);
for i = 0:9
    plot([sups(i+1),infs(i+1)],[i,i],'blue')
end

% plot samples
for i = 0:9
    scatter(ys(i+1,:), repmat(i,length(ys),1), '.k')
end

ylim([-1 9])
yticks(0:1:9)
yticklabels({'0','1','2','3','4','5','6','7','8','9'})
ylabel('Label')
xlabel('Prediction')
title('Set-based Prediction')

res = true;

% ------------------------------ END OF CODE ------------------------------
