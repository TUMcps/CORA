function res = test_nn_nnConvTranspose2DLayer(varargin)
% test_nn_nnConvTranspose2DLayer - tests constructor and evaluation of nnConvTranspose2DLayer
%
% Syntax:
%    res = test_nn_nnConvTranspose2DLayer()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% References:
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none

% See also: nnConvTranspose2DLayer

% Authors:       Benedikt Kellner
% Written:       01-December-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Validate parameters.
[options] = setDefaultValues({struct}, varargin);
options = nnHelper.validateNNoptions(options);

% --- 1. Test Constructor ---
% Standard constructor checks
layer = nnConvTranspose2DLayer([1 2 3; 4 5 6; 7 8 9]);
layer = nnConvTranspose2DLayer([1 2 3; 4 5 6; 7 8 9],1);
layer = nnConvTranspose2DLayer([1 2 3; 4 5 6; 7 8 9],1,[1 1 1 1]);
layer = nnConvTranspose2DLayer([1 2 3; 4 5 6; 7 8 9],1,[1 1 1 1],[2 2],[2 2]);
layer = nnConvTranspose2DLayer([1 2 3; 4 5 6; 7 8 9],1,[1 1 1 1],[2 2],[2 2],'testTransLayer');

% --- 2. Test Evaluate Numeric (Forward Pass) ---

% Parameters for a simple, verifiable transposed convolution
% Input image size: [2, 2, 1]
% Output channels (C_out): 2 (must match bias length)
% Filter size: [2, 2]
% Stride: [2, 2], Padding: [0 0 0 0]

% Define the 2 filters
W_filter1 = [1 2; 3 4]; 
W_filter2 = [5 6; 7 8]; 
b_transp = [1; 2]; % Bias for the 2 output channels (C_out = 2)

% Corrected W_transp definition to pass the nnConvTranspose2DLayer constructor check:
% W(k_h, k_w, C_out, C_in) -> [2, 2, 2, 1]
W_transp = zeros(2,2,2,1);
W_transp(:,:,1,1) = W_filter1; % Filter for 1st Output Channel
W_transp(:,:,2,1) = W_filter2; % Filter for 2nd Output Channel

layer = nnConvTranspose2DLayer(W_transp, b_transp, [0 0 0 0], [2 2]);
nn = neuralNetwork({layer});
n = 2; % Input image size: 2x2x1
nn.setInputSize([n,n,1]);

% Input vector x (2x2x1 = 4 neurons). Flattened column-major: [1 3; 2 4]
x = reshape([1 2; 3 4], [], 1); 

% Evaluate the layer
y = nn.evaluate(x,options);

% Expected Output Size: O = (I - 1) * S - P_total + F
% O_h = (2 - 1) * 2 - 0 + 2 = 4. Output size: 4x4x2 (32 neurons)

% Manual Calculation for Transposed Conv (Stride 2, Padding 0):
% y_out(i, j) = sum_k,l ( x(k,l) * W(i-k*S, j-l*S) ) + b
% Input x: [1 2; 3 4]. Filter 1: [1 2; 3 4], Bias 1: 1. Filter 2: [5 6; 7 8], Bias 2: 2.

% Output Channel 1: (Input values scaled by Filter 1 and scattered) + Bias 1
y1_map = [...
    1*1, 1*2, 2*1, 2*2;
    1*3, 1*4, 2*3, 2*4;
    3*1, 3*2, 4*1, 4*2;
    3*3, 3*4, 4*3, 4*4;
] + 1; 

% Output Channel 2: (Input values scaled by Filter 2 and scattered) + Bias 2
y2_map = [...
    1*5, 1*6, 2*5, 2*6;
    1*7, 1*8, 2*7, 2*8;
    3*5, 3*6, 4*5, 4*6;
    3*7, 3*8, 4*7, 4*8;
] + 2; 

y_true_map = cat(3, y1_map, y2_map);
y_true = reshape(y_true_map, [], 1); 

% Compare result
assert(all(abs(y - y_true) < 1e-3), 'Numeric evaluation failed for nnConvTranspose2DLayer.');

% --- 3. Check Interval Evaluation ---
% Create an interval object
inputInterval = interval(x - 0.1, x + 0.1);
% Evaluate
resInterval = nn.evaluate(inputInterval, options);
% Check if the resulting interval contains the numeric output (center)
assert(contains(resInterval, interval(y)), 'evaluateInterval failed: output interval does not contain the numeric result.');
% Check if the output is actually an interval object
assert(isa(resInterval, 'interval'), 'evaluateInterval failed: result is not an interval object.');


% --- 3. Check Zonotope Evaluation ---
% Use small perturbation for Zonotope test
perturbation = 0.001;
c1 = x - perturbation;
G1 = perturbation * eye(n*n);
Y1 = nn.evaluate(zonotope(c1,G1),options);
assert(contains(Y1,y), 'Zonotope bounds failed to contain center for nnConvTranspose2DLayer.');

c2 = x + 2*perturbation;
G2 = 2*perturbation * eye(n*n);
Y2 = nn.evaluate(zonotope(c2,G2),options);
assert(contains(Y2,y), 'Zonotope bounds failed to contain center for nnConvTranspose2DLayer.');

% --- 4. Check Zonotope Batch Evaluation 
nn.prepareForZonoBatchEval(x);
[cys,Gys] = nn.evaluateZonotopeBatch([c1 c2],cat(3,G1,G2));
assert(all(withinTol(cys(:,1),Y1.c,1e-10),'all'), 'Zonotope batch center 1 mismatch.');
assert(all(withinTol(cys(:,2),Y2.c,1e-10),'all'), 'Zonotope batch center 2 mismatch.');


assert(all(withinTol(Gys(:,:,1),Y1.G,1e-10),'all'), 'Zonotope batch generator 1 mismatch.');
assert(all(withinTol(Gys(:,:,2),Y2.G,1e-10),'all'), 'Zonotope batch generator 2 mismatch.');
assert(all(withinTol(Gys(:,:,1),Y1.G,1e-10),'all'), 'Zonotope batch generator 1 mismatch.');
assert(all(withinTol(Gys(:,:,2),Y2.G,1e-10),'all'), 'Zonotope batch generator 2 mismatch.');

% --- 6. Check PolyZonotope Evaluation ---
% Create a simple PolyZonotope
PZ = polyZonotope(zonotope(interval(x - 0.05, x + 0.05))); 
resPZ = nn.evaluate(PZ, options);
% Check containment
assert(contains(resPZ, y), 'evaluatePolyZonotope failed: result does not contain numeric input.');

% --- 7. Check Taylor Model Evaluation (taylm) ---
% Convert to Taylor Model
TM = taylm(zonotope(interval(x - 0.05, x + 0.05)));
resTM = nn.evaluate(TM, options);
% assert(in(interval(y), interval(resTM)), 'evaluateTaylm failed: result does not contain numeric input.');

% --- 8. Check Constrained Zonotope Evaluation (conZonotope) ---
CZ = conZonotope(zonotope(interval(x - 0.05, x + 0.05)));
resCZ = nn.evaluate(CZ, options);
% assert(contains(resCZ, y), 'evaluateConZonotope failed: result does not contain numeric input.');

% --- 9. Check evaluateSensitivity ---
% Create a dummy sensitivity matrix (Input Gradient). 
% Sensitivity is propagated from output to input (backwards), but evaluateSensitivity 
% assumes gradient w.r.t parameters or input?
% In nnLayer, evaluateSensitivity calls transconv2d. 
% For ConvTranspose2D, it should call conv2d (forward convolution).
% evaluateSensitivity is internal, so we check it via neuralNetwork.calcSensitivity.
% Output size of layer: 4x4x2 = 32.
outSize = layer.getOutputSize([n,n,1]);
numOut = prod(outSize);
dimIn = numel(x);

% Calculate sensitivity for input x
[S, ~] = nn.calcSensitivity(x, options);

% S should be (OutputDim x InputDim x BatchSize) = (32 x 4 x 1)
% Note: calcSensitivity initializes S as Identity(OutputDim), then backprops.
assert(size(S,1) == numOut, 'calcSensitivity failed: Output first dimension mismatch.');
assert(size(S,2) == dimIn, 'calcSensitivity failed: Output second dimension mismatch.');

% --- 10. Test Conversion to DLT ---
% Create a CORA Transposed Conv layer
% W input for CORA: (H, W, Out, In) -> (2, 2, 2, 1)
H = 2; W_dim = 2; Out = 2; In = 1;
W_in = rand(H, W_dim, Out, In, 'single');
b_in = rand(Out, 1, 'single');
padding_test = [0 0 0 0];
stride_test = [1 1];
dilation_test = [1 1];

cora_layer_conv = nnConvTranspose2DLayer(W_in, b_in, padding_test, stride_test, dilation_test, 'transp_conv');

% Create a minimal neuralNetwork containing this layer
nn_test = neuralNetwork({cora_layer_conv});
nn_test.neurons_in = [4 4 In]; % Dummy input size

% Convert to DLT
try
    dlt_net_test = nn_test.convertToDLToolboxNetwork();
    
    % Extract the layer
    if isa(dlt_net_test.Layers(1), 'nnet.cnn.layer.ImageInputLayer') || ...
       isa(dlt_net_test.Layers(1), 'nnet.cnn.layer.FeatureInputLayer')
       dlt_layer_test = dlt_net_test.Layers(2);
    else
       dlt_layer_test = dlt_net_test.Layers(1);
    end
    
    % Verify Class
    assert(isa(dlt_layer_test, 'nnet.cnn.layer.TransposedConvolution2DLayer'), 'Converted layer is not TransposedConv.');
    
    % Verify Weights
    % CORA Constructor takes (H, W, Out, In). Internal storage is (H, W, In, Out).
    % Converter reads Internal (H, W, In, Out), permutes back to (H, W, Out, In).
    % So DLT Weights should match the Input W_in exactly.
    assert(all(abs(dlt_layer_test.Weights - W_in) < 1e-6, 'all'), 'Converted weights mismatch.');
    
    % Verify Bias
    assert(all(abs(reshape(dlt_layer_test.Bias, [], 1) - b_in) < 1e-6, 'all'), 'Converted bias mismatch.');
catch ME
    CORAwarning('CORA:nn',sprintf('DLT Conversion Test Failed or DLT not available: %s', ME.message));
    % Rethrow if it's an assertion error
    if strcmp(ME.identifier, 'MATLAB:assertion:failed')
        rethrow(ME);
    end
end

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
