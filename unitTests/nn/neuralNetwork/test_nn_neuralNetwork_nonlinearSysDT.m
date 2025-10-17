function res = test_nn_neuralNetwork_nonlinearSysDT
% test_nn_neuralNetwork_nonlinearSysDT - unit test for transforming a  
%       neural network object to a nonlinearSysDT object
%
% Syntax:
%    res = test_nn_neuralNetwork_nonlinearSysDT
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% References:
%    [1] Laura Luetzow, Michael Eichelbeck, Mykel Kochenderfer, and
%        Matthias Althoff. "Zono-Conformal Prediction: Zonotope-Based
%        Uncertainty Quantification for Regression and Classification
%        Tasks," arXiv, 2025.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: setPredictor

% Authors:       Laura Luetzow
% Written:       08-October-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

for i = 1:3
    % single linear layer
    dim_x = randi([1 5],1);
    dim_y = randi([1 5],1);
    W = randn(dim_y,dim_x);
    nn = neuralNetwork({nnLinearLayer(W)});
    sys = nonlinearSysDT(nn);

    % check dimensions
    assert(sys.nrOfDims == dim_x)
    assert(sys.nrOfInputs == dim_y) % all biases are uncertain
    assert(sys.nrOfOutputs == dim_y)
    x_i = randn(dim_x,1);
    u_i = randn(dim_y,1);

    % check dynamics
    assert(all(withinTol(sys.mFile(x_i,u_i), x_i, 1e-6),'all'))
    assert(all(withinTol(sys.out_mFile(x_i,u_i), W*x_i + u_i, 1e-6),'all'))
    assert(all(withinTol(sys.jacobian(x_i,u_i), eye(dim_x), 1e-6),'all'))
    assert(all(withinTol(sys.out_jacobian(x_i,u_i), W, 1e-6),'all'))

    % add nonlinear layer
    nn = neuralNetwork({nnLinearLayer(W); nnReLULayer});
    sys = nonlinearSysDT(nn);

    % check dimensions
    assert(sys.nrOfDims == dim_x)
    assert(sys.nrOfInputs == dim_y) % all biases are uncertain
    assert(sys.nrOfOutputs == dim_y)
    x_i = randn(dim_x,1);
    u_i = randn(dim_y,1);

    % check dynamics
    assert(all(withinTol(sys.mFile(x_i,u_i), x_i, 1e-6),'all'))
    y_i = W*x_i + u_i;
    y_i(y_i < 0) = 0;
    assert(all(withinTol(sys.out_mFile(x_i,u_i), y_i, 1e-6),'all'))
    assert(all(withinTol(sys.jacobian(x_i,u_i), eye(dim_x), 1e-6),'all'))

    % jacobian is computed setting u_i = 0!!!
    J_i = W;
    J_i(W*x_i < 0,:) = 0;
    assert(all(withinTol(sys.out_jacobian(x_i,u_i), J_i, 1e-6),'all'))

    % add two more layers
    dim_y2 = randi([1 5],1);
    W2 = randn(dim_y2,dim_y);
    nn = neuralNetwork({nnLinearLayer(W); nnReLULayer(); ...
        nnLinearLayer(W2); nnReLULayer()});
    sys = nonlinearSysDT(nn);

    % check dimensions
    assert(sys.nrOfDims == dim_x)
    assert(sys.nrOfInputs == dim_y+dim_y2) % all biases are uncertain
    assert(sys.nrOfOutputs == dim_y2)
    x_i = randn(dim_x,1);
    u_i = randn(dim_y+dim_y2,1);

    % check dynamics
    assert(all(withinTol(sys.mFile(x_i,u_i), x_i, 1e-6),'all'))
    y_i = W*x_i + u_i(1:dim_y);
    y_i(y_i < 0) = 0;
    y_i2 = W2*y_i + u_i(dim_y+1:end);
    y_i2(y_i2 < 0) = 0;
    assert(all(withinTol(sys.out_mFile(x_i,u_i), y_i2, 1e-6),'all'))
    assert(all(withinTol(sys.jacobian(x_i,u_i), eye(dim_x), 1e-6),'all'))
    J_i = W;
    J_i(W*x_i < 0,:) = 0;
    J_i2 = W2*J_i;
    J_i2(W2*J_i*x_i < 0,:) = 0;
    assert(all(withinTol(sys.out_jacobian(x_i,u_i), J_i2, 1e-6),'all'))
end

res = true;
end

% ------------------------------ END OF CODE ------------------------------
