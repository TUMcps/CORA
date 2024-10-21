function res = test_nn_nnSGDOptimizer_step()
% test_nn_nnSGDOptimizer_step - unit test function for 
%     nnSGDOptimizer.step: compare update step with MATLAB Deep
%     Learning Toolbox
%
% Syntax:
%    res = test_nn_nnSGDOptimizer_step()
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
% See also: nnSGDOptimizer

% Authors:       Lukas Koller
% Written:       10-March-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
% assume true
res = true;

rng('default')

% Specify size of linear layer.
n_in = 1;
n_out = 1;

% Specify number of iterations.
N = 10;

% Construct a linear layer with random weights.
W = 2*rand([n_out n_in]) - 1; % Sample random weight from [-1,1]
linlayer = nnLinearLayer(W,0,'linlayer');
nn = neuralNetwork({linlayer});

% Sample random gradients.
dW = cell(N,1);
for i=1:N
    dW{i} = 2*rand([n_out n_in]) - 1;
end

% Sample random initial gradient velocity.
vel = 2*rand([n_out n_in]) - 1;

% Specify time step.
timestep = 42;

% Specify SGD parameters.
lr = 1e-1;
momentum = 0.9;

% CORA IMPLEMENTATION

% Get default options.
opts.nn = struct('poly_method','bounds');
opts = nnHelper.validateNNoptions(opts,true);

% Construct optimizer.
optim = nnSGDOptimizer(lr,momentum);
% Init optimizer.
optim.deleteGrad(nn,opts);
optim.timestep = timestep;

nextWCora = cell(N+1,1);
velCora = cell(N+1,1);
% Set initial weights.
nextWCora{1} = W;
% Set initial velocity.
velCora{1} = vel;
for i=1:N
    % Set gradient.
    linlayer.backprop.grad.('W') = dW{i};
    % Set velocity.
    linlayer.backprop.vel.('W') = velCora{i};

    % Update parameters.
    optim.step(nn,opts);
    
    % Extract updated weights.
    nextWCora{i+1} = linlayer.W;
    % Extract updated gradient velocity.
    velCora{i+1} = linlayer.backprop.vel.('W');
end

% DEEP LEARNING TOOLBOX IMPLEMENTATION

nextWDlt = cell(N+1,1);
velDlt = cell(N+1,1);
% Set initial weights.
nextWDlt{1} = W;
% Set initial velocity.
velDlt{1} = vel;
for i=1:N
    [nextWDlt{i+1},velDlt{i+1}] = sgdmupdate(nextWDlt{i},dW{i},velDlt{i},lr,momentum);
end

% COMPARISON
for i=1:N+1
    assertLoop(aux_eq(nextWCora{i},nextWDlt{i}),i);
    assertLoop(aux_eq(velCora{i},velDlt{i}),i);
end

end


% Auxiliary functions -----------------------------------------------------

function r = aux_eq(a,b)
    % r = all(abs(a - b) <= eps('like',a),'all');
    r = all(abs(a - b) <= 1e-10,'all');
end

% ------------------------------ END OF CODE ------------------------------
