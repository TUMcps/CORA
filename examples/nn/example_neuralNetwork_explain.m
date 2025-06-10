function example_neuralNetwork_explain(XTest, YTest)
% example_neuralNetwork_explain - example for xai in cora
%
% Syntax:
%    res = example_neuralNetwork_explain()
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
% Written:       27-May-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

seed = randi(1000)
rng(seed)

epsilon = 0.01;

% Load model and data -----------------------------------------------------

disp("Reading model and data...") % (onnx requires us to specify input format)
modelfile = "mnist_sigmoid_6_200.onnx";
% modelfile = "mnist_relu_6_200.onnx";
% modelfile = "convMedGRELU__Point.onnx";
nn = neuralNetwork.readONNXNetwork(modelfile, false, 'BCSS');
% modelfile = "fnn_mnist_3x200.onnx";
% nn = neuralNetwork.readONNXNetwork(modelfile, true, 'BC');
% modelfile = "mnist-net_256x4.onnx";
% nn = neuralNetwork.readONNXNetwork(modelfile, false, 'CSS');

% reformulate network to obtain: linear -> activation -> linear -> ...
nn = nn.getNormalForm();

% load data
if nargin < 2
    [XTest,YTest, ~] = digitTest4DArrayData;
end

disp("Done.")

N = 1; % do it for one image

% Search for correctly predicted images -----------------------------------

images = aux_searchCorrectlyPredictedImages(N, nn, XTest, YTest);

% Explain each image ------------------------------------------------------

for i=images
    fprintf("Explain image: %i\n", i);
    
    % reset network before verification
    nn.reset();

    % get image i
    x = XTest(:, :, 1, i);
    target = double(YTest(i));

    % reshape to column vector
    x = reshape(x, [], 1);
    x = double(x);

    % call explain
    timerVal = tic;
    nn.explain(x,target,epsilon,'Method','abstract+refine','Verbose',true,'InputSize',[28,28,1]);
    toc(timerVal);
   
end

end


% Auxiliary functions -----------------------------------------------------

function images = aux_searchCorrectlyPredictedImages(N, nn, XTest, YTest)

images = [];
while(length(images) ~= N)

    % select random image
    i = randi(length(YTest));
    x = XTest(:, :, 1, i);
    x = double(x);
    target = YTest(i);
    fprintf("Label: %i - ", target)

    % reshape to column vector
    x = reshape(x, [], 1);

    % evaluate output
    y = nn.evaluate(x);
    [~, pred] = max(y);
    
    % check with target
    correct = pred == double(target);
    if correct
        images = [images, i];
        fprintf("added.\n")
    else 
        fprintf("misclassified image.\n")
    end
end
disp(" ")

end


% ------------------------------ END OF CODE ------------------------------
