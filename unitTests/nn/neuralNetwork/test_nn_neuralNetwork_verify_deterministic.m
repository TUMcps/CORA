function res = test_nn_neuralNetwork_verify_deterministic()
% test_nn_neuralNetwork_verify_deterministic - unit test for
%    neuralNetwork/verify using hand-constructed networks with known results
%
% Syntax:
%    res = test_nn_neuralNetwork_verify_deterministic()
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
% See also: neuralNetwork/verify

% Authors:       Benedikt Kellner, Lukas Koller
% Written:       12-March-2026
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% Set seed for reproducibility.
rng(1);

% Timeout for all tests (small networks should be instant).
timeout = 5;
verbose = false;

% --- Build test networks --------------------------------------------------

% Identity network: y = x (1-in/1-out).
nnId = neuralNetwork({nnLinearLayer(1, 0)});

% ReLU network: y = ReLU(x1) + ReLU(x2) (2-in/1-out).
nnReLU = neuralNetwork({ ...
    nnLinearLayer(eye(2), [0;0]); ...
    nnReLULayer(); ...
    nnLinearLayer([1 1], 0) ...
});

% --- Helper to build baseline options -------------------------------------

options = aux_baselineOptions();

% =========================================================================
% Test 1: Identity network, trivially VERIFIED
% =========================================================================
x = 0; r = 1; A = 1; b = 2; safeSet = true;

% 1a: naive refinement
opts = options;
opts.nn.refinement_method = 'naive';
[verifRes,x_,y_] = nnId.verify(x,r,A,b,safeSet,opts,timeout,verbose);
assert(strcmp(verifRes.str,'VERIFIED'));
assert(isempty(x_) & isempty(y_));

% 1b: zonotack refinement
opts.nn.refinement_method = 'zonotack';
[verifRes,x_,y_] = nnId.verify(x,r,A,b,safeSet,opts,timeout,verbose);
assert(strcmp(verifRes.str,'VERIFIED'));
assert(isempty(x_) & isempty(y_));

% =========================================================================
% Test 2: Identity network, trivially COUNTEREXAMPLE
% =========================================================================
x = 0; r = 1; A = 1; b = 0; safeSet = true;
[verifRes,x_,y_] = nnId.verify(x,r,A,b,safeSet,options,timeout,verbose);
assert(strcmp(verifRes.str,'COUNTEREXAMPLE'));
aux_validateCounterexample(nnId,x,r,A,b,safeSet,x_,y_);

% =========================================================================
% Test 3: ReLU network, VERIFIED
% =========================================================================
x = [0;0]; r = [0.5;0.5]; A = 1; b = 3; safeSet = true;
[verifRes,x_,y_] = nnReLU.verify(x,r,A,b,safeSet,options,timeout,verbose);
assert(strcmp(verifRes.str,'VERIFIED'));
assert(isempty(x_) & isempty(y_));

% =========================================================================
% Test 4: ReLU network, COUNTEREXAMPLE
% =========================================================================
x = [0;0]; r = [0.5;0.5]; A = 1; b = 0.3; safeSet = true;
[verifRes,x_,y_] = nnReLU.verify(x,r,A,b,safeSet,options,timeout,verbose);
assert(strcmp(verifRes.str,'COUNTEREXAMPLE'));
aux_validateCounterexample(nnReLU,x,r,A,b,safeSet,x_,y_);

% =========================================================================
% Test 5: Unsafe set specification
% =========================================================================
x = 0; r = 1; A = 1; b = 0.5; safeSet = false;
% Unsafe if y <= 0.5; output in [-1,1], so counterexample exists.
[verifRes,x_,y_] = nnId.verify(x,r,A,b,safeSet,options,timeout,verbose);
assert(strcmp(verifRes.str,'COUNTEREXAMPLE'));
aux_validateCounterexample(nnId,x,r,A,b,safeSet,x_,y_);

% =========================================================================
% Test 6: All falsification methods on known-COUNTEREXAMPLE case
% =========================================================================
x = 0; r = 1; A = 1; b = 0; safeSet = true;
falsifMethods = {'center','fgsm','zonotack'};
for i = 1:length(falsifMethods)
    opts = options;
    opts.nn.falsification_method = falsifMethods{i};
    [verifRes,x_,y_] = nnId.verify(x,r,A,b,safeSet,opts,timeout,verbose);
    assert(strcmp(verifRes.str,'COUNTEREXAMPLE'), ...
        sprintf('Falsification method %s failed.',falsifMethods{i}));
    aux_validateCounterexample(nnId,x,r,A,b,safeSet,x_,y_);
end

% =========================================================================
% Test 7: All refinement methods on known-VERIFIED case
% =========================================================================
x = 0; r = 1; A = 1; b = 2; safeSet = true;
refMethods = {'naive','zonotack','zonotack-layerwise'};
for i = 1:length(refMethods)
    opts = options;
    opts.nn.refinement_method = refMethods{i};
    [verifRes,x_,y_] = nnId.verify(x,r,A,b,safeSet,opts,timeout,verbose);
    assert(strcmp(verifRes.str,'VERIFIED'), ...
        sprintf('Refinement method %s failed.',refMethods{i}));
    assertLoop(isempty(x_) & isempty(y_),i);
end

% =========================================================================
% Test 8: Queue strategies on known-VERIFIED case
% =========================================================================
x = 0; r = 1; A = 1; b = 2; safeSet = true;
deqTypes = {'front','half-half'};
enqTypes = {'append','prepend'};
for i = 1:length(deqTypes)
    for j = 1:length(enqTypes)
        opts = options;
        opts.nn.verify_dequeue_type = deqTypes{i};
        opts.nn.verify_enqueue_type = enqTypes{j};
        [verifRes,~,~] = nnId.verify(x,r,A,b,safeSet,opts,timeout,verbose);
        assertLoop(strcmp(verifRes.str,'VERIFIED'),i,j);
    end
end

% =========================================================================
% Test 9: Splitting parameters on known-VERIFIED ReLU case
% =========================================================================
x = [0;0]; r = [0.5;0.5]; A = 1; b = 3; safeSet = true;

splitConfigs = {
    struct('num_splits',1,'num_dimensions',0,'num_neuron_splits',0, ...
           'input_xor_neuron_splitting',false,'add_orth_neuron_splits',false);
    struct('num_splits',3,'num_dimensions',2,'num_neuron_splits',1, ...
           'input_xor_neuron_splitting',false,'add_orth_neuron_splits',false);
    struct('num_splits',2,'num_dimensions',1,'num_neuron_splits',0, ...
           'input_xor_neuron_splitting',true,'add_orth_neuron_splits',false);
    struct('num_splits',2,'num_dimensions',1,'num_neuron_splits',1, ...
           'input_xor_neuron_splitting',false,'add_orth_neuron_splits',true);
};
for i = 1:length(splitConfigs)
    opts = options;
    flds = fieldnames(splitConfigs{i});
    for f = 1:length(flds)
        opts.nn.(flds{f}) = splitConfigs{i}.(flds{f});
    end
    [verifRes,~,~] = nnReLU.verify(x,r,A,b,safeSet,opts,timeout,verbose);
    assertLoop(~strcmp(verifRes.str,'COUNTEREXAMPLE'),i);
end

% =========================================================================
% Test 10: Heuristics
% =========================================================================

% 10a: input_generator_heuristic
x = [0;0]; r = [0.5;0.5]; A = 1; b = 3; safeSet = true;
inGenHeurs = {'most-sensitive-input-radius','zono-norm-gradient'};
for i = 1:length(inGenHeurs)
    opts = options;
    opts.nn.input_generator_heuristic = inGenHeurs{i};
    [verifRes,~,~] = nnReLU.verify(x,r,A,b,safeSet,opts,timeout,verbose);
    assertLoop(~strcmp(verifRes.str,'COUNTEREXAMPLE'),i);
end

% 10b: input_split_heuristic
inSplitHeurs = {'most-sensitive-input-radius','zono-norm-gradient'};
for i = 1:length(inSplitHeurs)
    opts = options;
    opts.nn.input_split_heuristic = inSplitHeurs{i};
    [verifRes,~,~] = nnReLU.verify(x,r,A,b,safeSet,opts,timeout,verbose);
    assertLoop(~strcmp(verifRes.str,'COUNTEREXAMPLE'),i);
end

% 10c: neuron_split_heuristic (requires num_neuron_splits >= 1)
neuronHeurs = {'least-unstable','least-unstable-gradient', ...
    'most-sensitive-approx-error','most-sensitive-input-radius', ...
    'zono-norm-gradient'};
for i = 1:length(neuronHeurs)
    opts = options;
    opts.nn.num_neuron_splits = 1;
    opts.nn.neuron_split_heuristic = neuronHeurs{i};
    [verifRes,~,~] = nnReLU.verify(x,r,A,b,safeSet,opts,timeout,verbose);
    assertLoop(~strcmp(verifRes.str,'COUNTEREXAMPLE'),i);
end

% 10d: relu_constraint_heuristic (requires num_relu_constraints >= 1)
reluHeurs = {'least-unstable','least-unstable-gradient', ...
    'most-sensitive-approx-error','most-sensitive-input-radius', ...
    'zono-norm-gradient'};
for i = 1:length(reluHeurs)
    opts = options;
    opts.nn.num_relu_constraints = 3;
    opts.nn.relu_constraint_heuristic = reluHeurs{i};
    [verifRes,~,~] = nnReLU.verify(x,r,A,b,safeSet,opts,timeout,verbose);
    assertLoop(~strcmp(verifRes.str,'COUNTEREXAMPLE'),i);
end

% =========================================================================
% Test 11: Approximation error options
% =========================================================================
x = [0;0]; r = [0.5;0.5]; A = 1; b = 3; safeSet = true;

% 11a: approx_error_order
for i = 1:2
    orders = {'length','sensitivity*length'};
    opts = options;
    opts.nn.approx_error_order = orders{i};
    [verifRes,~,~] = nnReLU.verify(x,r,A,b,safeSet,opts,timeout,verbose);
    assertLoop(~strcmp(verifRes.str,'COUNTEREXAMPLE'),i);
end

% 11b: train.num_approx_err
for i = 1:2
    vals = {inf, 5};
    opts = options;
    opts.nn.train.num_approx_err = vals{i};
    [verifRes,~,~] = nnReLU.verify(x,r,A,b,safeSet,opts,timeout,verbose);
    assertLoop(~strcmp(verifRes.str,'COUNTEREXAMPLE'),i);
end

% 11c: train.num_init_gens
for i = 1:2
    vals = {inf, 2};
    opts = options;
    opts.nn.train.num_init_gens = vals{i};
    [verifRes,~,~] = nnReLU.verify(x,r,A,b,safeSet,opts,timeout,verbose);
    assertLoop(~strcmp(verifRes.str,'COUNTEREXAMPLE'),i);
end

% =========================================================================
% Test 12: Interval center mode
% =========================================================================
x = [0;0]; r = [0.5;0.5]; A = 1; b = 3; safeSet = true;
for i = 1:2
    vals = {true, false};
    opts = options;
    opts.nn.interval_center = vals{i};
    [verifRes,~,~] = nnReLU.verify(x,r,A,b,safeSet,opts,timeout,verbose);
    assertLoop(~strcmp(verifRes.str,'COUNTEREXAMPLE'),i);
end

% =========================================================================
% Test 13: ReLU constraints
% =========================================================================
x = [0;0]; r = [0.5;0.5]; A = 1; b = 3; safeSet = true;
reluConsts = [0, 3, inf];
for i = 1:length(reluConsts)
    opts = options;
    opts.nn.num_relu_constraints = reluConsts(i);
    [verifRes,~,~] = nnReLU.verify(x,r,A,b,safeSet,opts,timeout,verbose);
    assertLoop(~strcmp(verifRes.str,'COUNTEREXAMPLE'),i);
end

% =========================================================================
% Test 14: max_verif_iter = 1
% =========================================================================
x = [0;0]; r = [0.5;0.5]; A = 1; b = 3; safeSet = true;
opts = options;
opts.nn.max_verif_iter = 1;
[verifRes,~,~] = nnReLU.verify(x,r,A,b,safeSet,opts,timeout,verbose);
% With max_verif_iter=1, result may be VERIFIED or UNKNOWN, never COUNTEREXAMPLE.
assert(~strcmp(verifRes.str,'COUNTEREXAMPLE'));

end


% Auxiliary functions -----------------------------------------------------

function options = aux_baselineOptions()
    options.nn = struct(...
        'use_approx_error',true,...
        'poly_method','bounds',...
        'train',struct(...
            'backprop',false,...
            'mini_batch_size',2^8 ...
        ) ...
    );
    options = nnHelper.validateNNoptions(options,true);
    options.nn.interval_center = false;
    % Use zonotack for both falsification and refinement by default.
    options.nn.falsification_method = 'zonotack';
    options.nn.refinement_method = 'zonotack';
    options.nn.num_splits = 2;
    options.nn.num_dimensions = 1;
    options.nn.num_neuron_splits = 0;
end

function aux_validateCounterexample(nn,x,r,A,b,safeSet,x_,y_)
    assert(~isempty(x_) & ~isempty(y_));
    % Check input bounds.
    tol = 1e-6;
    assert(all(x_ >= x - r - tol,'all') & all(x_ <= x + r + tol,'all'), ...
        'Counterexample x_ out of input bounds.');
    % Check output matches network evaluation.
    yi = nn.evaluate(x_);
    assert(all(abs(y_ - yi) <= 1e-7,'all'), ...
        'Counterexample y_ does not match nn.evaluate(x_).');
    % Check specification is violated.
    if safeSet
        assert(any(A*yi >= b,1), ...
            'Counterexample does not violate safe-set specification.');
    else
        assert(all(A*yi <= b,1), ...
            'Counterexample does not violate unsafe-set specification.');
    end
end

% ------------------------------ END OF CODE ------------------------------
