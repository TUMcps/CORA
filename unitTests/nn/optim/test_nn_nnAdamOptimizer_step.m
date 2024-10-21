function res = test_nn_nnAdamOptimizer_step()
% test_nn_nnAdamOptimizer_step - unit test function for 
%     nnAdamOptimizer.step: compare update step with MATLAB Deep
%     Learning Toolbox
%
% Syntax:
%    res = test_nn_nnAdamOptimizer_step()
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
% See also: nnAdamOptimizer

% Authors:       Lukas Koller
% Written:       10-March-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

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

% Sample random inital moment estimates.
avgGrad = 2*rand([n_out n_in]) - 1;
avgSqGrad = 2*rand([n_out n_in]) - 1;

% Specify time step.
timestep = 42;

% Specify Adam parameters.
lr = 1e-1;
beta1 = 0.9;
beta2 = 0.999;
epsilon = 1e-8;

% CORA IMPLEMENTATION

% Get default options.
opts.nn = struct('poly_method','bounds');
opts = nnHelper.validateNNoptions(opts,true);

% Construct optimizer.
optim = nnAdamOptimizer(lr,beta1,beta2,epsilon);
% Init optimizer.
optim.deleteGrad(nn,opts);
optim.timestep = timestep;

nextWCora = cell(N+1,1);
avgGradCora = cell(N+1,1);
avgSqGradCora = cell(N+1,1);
% Set initial weights.
nextWCora{1} = W;
% Set initial moment estimates.
avgGradCora{1} = avgGrad;
avgSqGradCora{1} = avgSqGrad;
for i=1:N
    % Set gradient.
    linlayer.backprop.grad.('W') = dW{i};
    % Set moment estimates.
    linlayer.backprop.mt.('W') = avgGradCora{i};
    linlayer.backprop.vt.('W') = avgSqGradCora{i};

    % Update parameters.
    optim = optim.step(nn,opts);
    
    % Extract updated weights.
    nextWCora{i+1} = linlayer.W;
    % Extract updated moment estimates.
    avgGradCora{i+1} = linlayer.backprop.mt.('W');
    avgSqGradCora{i+1} = linlayer.backprop.vt.('W');
end

% DEEP LEARNING TOOLBOX IMPLEMENTATION

nextWDlt = cell(N+1,1);
avgGradDlt = cell(N+1,1);
avgSqGradDlt = cell(N+1,1);
% Set initial weights.
nextWDlt{1} = W;
% Set initial moment estimates.
avgGradDlt{1} = avgGrad;
avgSqGradDlt{1} = avgSqGrad;
for i=1:N
    % There's a bug in the Deep Learning Toolbox: the bias correction is
    % not correct applied and does not consider the epsilon. Therefore, we
    % manually adjust the epsilon to compensate for the bug.
    epsilonBc = epsilon.*sqrt(1 - beta2.^(timestep+i));
    [nextWDlt{i+1},avgGradDlt{i+1},avgSqGradDlt{i+1}] = ...
        adamupdate(nextWDlt{i},dW{i},avgGradDlt{i},avgSqGradDlt{i}, ...
            timestep+i,lr,beta1,beta2,epsilonBc);
end

% COMPARISON
for i=1:N+1
    assertLoop(aux_eq(nextWCora{i},nextWDlt{i}),i);
    assertLoop(aux_eq(avgGradCora{i},avgGradDlt{i}),i);
    assertLoop(aux_eq(avgSqGradCora{i},avgSqGradDlt{i}),i);
end

% test completed
res = true;

end


% Auxiliary functions -----------------------------------------------------

function r = aux_eq(a,b)
    % r = all(abs(a - b) <= eps('like',a),'all');
    r = all(abs(a - b) <= 1e-10,'all');
end

% ------------------------------ END OF CODE ------------------------------
