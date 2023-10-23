function [res] = test_nn_neuralNetwork_getNormalForm()
% test_nn_neuralNetwork_getNormalForm - tests the
%    neuralNetwork/getNormalForm function
%
% Syntax:
%    res = test_nn_neuralNetwork_getNormalForm()
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
% Written:       14-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% test small network 
nn = neuralNetwork({ ...
    nnElementwiseAffineLayer(4, 2); ...
    nnLinearLayer(rand(3, 2), rand(3, 1)); ...
    nnSigmoidLayer(); ...
});
nn_normal = nn.getNormalForm();
res = res & length(nn_normal.layers) == 2;
res = res & aux_test_equality(nn, nn_normal);

nn = neuralNetwork({ ...
    nnElementwiseAffineLayer(1, -4); ...
    nnElementwiseAffineLayer(2, 0); ...
    nnLinearLayer(rand(3, 2), rand(3, 1)); ...
    nnSigmoidLayer(); ...
});
nn_normal = nn.getNormalForm();
res = res & length(nn_normal.layers) == 2;
res = res & aux_test_equality(nn, nn_normal);

% check non-scalar elementwise operations
nn = neuralNetwork({ ...
    nnLinearLayer(rand(2, 2), rand(2, 1));
    nnElementwiseAffineLayer([2;2], -4); ...
    nnElementwiseAffineLayer(2, 0); ...
    nnLinearLayer(rand(3, 2), rand(3, 1)); ...
    nnElementwiseAffineLayer(4, [1;2;3]); ...
    nnElementwiseAffineLayer([2;3;1], [1;2;3]); ...
    nnLinearLayer(rand(3, 3), rand(3, 1)); ...
    nnSigmoidLayer(); ...
});
nn_normal = nn.getNormalForm();
res = res & length(nn_normal.layers) == 2;
res = res & aux_test_equality(nn, nn_normal);

% test large network
nn = neuralNetwork({ ...
    nnElementwiseAffineLayer(4, 2); ...
    nnLinearLayer(rand(3, 2), rand(3, 1)); ...
    nnSigmoidLayer(); ...
    nnLinearLayer(rand(4, 3), rand(4, 1));
    nnTanhLayer(); ...
    nnIdentityLayer();
    nnLinearLayer(rand(2, 4), rand(2, 1));
    nnElementwiseAffineLayer(0.1, 2);
    nnReLULayer(); ...
});
nn_normal = nn.getNormalForm();
res = res & length(nn_normal.layers) == 6;
res = res & aux_test_equality(nn, nn_normal);


% test empty network
nn = neuralNetwork({});
nn_normal = nn.getNormalForm();
res = res & isempty(nn_normal.layers);
res = res & aux_test_equality(nn, nn_normal);


end


% Auxiliary functions -----------------------------------------------------

function res = aux_test_equality(nn, nn_normal)
    n = nn.neurons_in;
    N = 10;
    xs = rand(n, N);

    ys = nn.evaluate(xs);
    ys_normal = nn_normal.evaluate(xs);

    res = all(withinTol(ys, ys_normal), 'all');
end

% ------------------------------ END OF CODE ------------------------------
