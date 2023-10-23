function res = test_nn_zonotope_relu()
% test_nn_zonotope_relu - tests nn with relu activation using zonotopes
%
% Syntax:
%    res = test_nn_zonotope_relu()
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

% Authors:       Tobias Ladner
% Written:       24-June-2022
% Last update:   30-November-2022 (point inclusion check)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% load W, b, input_ref, output_ref (from previous model)
load('model_test_nn_zonotope_relu.mat');

% build the new layer based model
layers = {};
for i = 1:length(W)
    % concat layers
    layers = [; ...
        layers; ...
        {nnLinearLayer(W{i}, b{i})}; ...
        {nnReLULayer()}; ...
        ];
end
nn_new = neuralNetwork(layers);

evParams = struct;
evParams.poly_method = 'singh';

% calculate output
output_new = nn_new.evaluate(input_ref, evParams);
res = isequal(output_ref, output_new, 1e-15);

% compute points
points = [input_ref.randPoint(500), input_ref.randPoint(500, 'extreme')];
points = nn_new.evaluate(points);
res = res && all(contains(output_ref, points));

end

% ------------------------------ END OF CODE ------------------------------
