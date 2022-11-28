function res = test_nn_taylm_tanh()
% test_nn_taylm_tanh - tests nn with tanh activation using taylor models 
%
% Syntax:  
%    res = test_nn_taylm_tanh()
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

% Author:       Tobias Ladner
% Written:      17-November-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% load W, b, input_ref, output_ref (from previous model)
model = "model_test_nn_taym_tanh.mat";
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
    is_close(infimum(interval(output_ref)), infimum(interval(output_new))); ...
    is_close(supremum(interval(output_ref)), supremum(interval(output_new))); ...
    ];
res = all(res);
if ~res
    disp(["Failed!", model])
end

end

%------------- END OF CODE --------------