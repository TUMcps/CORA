function res = testFlaky_nn_polyZonotope_cnn_random()
% testFlaky_nn_polyZonotope_cnn_random - test a random CNN
%
%
% Syntax:
%    res = testFlaky_nn_polyZonotope_cnn_random
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
% Written:       17-August-2022
% Last update:   30-January-2024 (TL, flaky until merge with nn-gd)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% INPUT SET

img_h = 28;
img_w = 28;
c_in = 3;
c_hidden = 4;
c_out = 2;
S_in = polyZonotope.generateRandom('Dimension',img_h*img_w*c_in,'NrGenerators',15);

% CREATE NETWORK

% random convolutional weights
w = [3, 3]; % window size (h, w)
CW1 = rand(w(1), w(2), c_in, c_hidden) * 2 -1;
Cb1 = rand(c_hidden, 1) * 2 -1;
CW2 = rand(w(1), w(2), c_hidden, c_out) * 2 -1;
Cb2 = rand(c_out, 1) * 2 -1;

layers = {
    nnConv2DLayer(CW1, Cb1);
    nnAvgPool2DLayer([2, 2]);
    nnConv2DLayer(CW2 , Cb2);
    nnMaxPool2DLayer([2, 2]);
};
nn_cora = neuralNetwork(layers);

verbose = true;
nn_cora.setInputSize([img_h,img_w,c_in], verbose);

% RUN EVALUATE

evParams = struct();
evParams.bound_approx = true;
evParams.num_generators = 100;

S_out = nn_cora.evaluate(S_in, evParams);

% test sensitivity
nn_cora.calcSensitivity(S_in.c);

% TEST FOR POINTS

P_in = [S_in.randPoint(100), S_in.randPoint(50, 'extreme')];
P_out = nn_cora.evaluate(P_in);
res = all(zonotope(S_out).contains(P_out));

end

% ------------------------------ END OF CODE ------------------------------
