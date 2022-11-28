function res = test_nn_conZonotope_ReLU_exact()
% test_nn_conZonotope_ReLU_exact - tests nn with relu activation using
% conZonotopes (exact)
%
% Syntax:  
%    res = test_nn_conZonotope_ReLU_exact()
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
% Written:      24-June-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% load W, b, input_ref, output_ref (from previous model)
model = "model_test_nn_conZonotope_ReLU_exact.mat";
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
    is_close(output_ref.Z, output_new.Z); ...
    is_close(output_ref.A, output_new.A); ...
    is_close(output_ref.b, output_new.b); ...
    ];
res = all(res);
if ~res
    disp(["Failed!", model])
end

end

%------------- END OF CODE --------------