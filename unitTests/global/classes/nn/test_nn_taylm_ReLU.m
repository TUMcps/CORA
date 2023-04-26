function res = test_nn_taylm_ReLU()
% test_nn_taylm_ReLU - tests nn with relu activation using taylor models 
%
% Syntax:  
%    res = test_nn_taylm_ReLU()
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

% Author:       Tobias Ladner
% Written:      24-June-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% load W, b, input_ref, output_ref (from previous model)
model = "model_test_nn_taym_ReLU.mat";
load(model)

% build the new layer based model
layers = {};
for i = 1:length(W)
    % concat layers
    activation_layer = nnActivationLayer.instantiateFromString(evParams.activation);
    layers = [; ...
        layers; ...
        {nnLinearLayer(W{i}, b{i})}; ...
        {activation_layer}; ...
        ];
end
nn_new = neuralNetwork(layers);

% calculate output
output_new = nn_new.evaluate(input_ref, evParams);

res = [; ...
    is_close(interval(output_ref).inf, interval(output_new).inf); ...
    is_close(interval(output_ref).sup, interval(output_new).sup);
    ];
res = all(res);
if ~res
    disp(["Failed!", model])
end

end

%------------- END OF CODE --------------