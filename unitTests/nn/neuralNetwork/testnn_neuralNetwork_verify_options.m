function res = testnn_neuralNetwork_verify_options()
% testnn_neuralNetwork_verify_options - systematic test of all algorithm
%    options for neuralNetwork/verify using a real ONNX model
%
% Syntax:
%    res = testnn_neuralNetwork_verify_options()
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

% Toggle verbose verification output.
verbose = false;

% Timeout per verification call.
timeout = 1;

% --- Load known-UNSAT model (prop_1) ------------------------------------
model1Path = [CORAROOT '/models/Cora/nn/ACASXU_run2a_1_2_batch_2000.onnx'];
prop1Path = [CORAROOT '/models/Cora/nn/prop_1.vnnlib'];
[nn1,opts1,x1,r1,A1,b1,safeSet1] = aux_readNetworkAndOptions( ...
    model1Path,prop1Path);

% --- Load known-SAT model (prop_2) --------------------------------------
prop2Path = [CORAROOT '/models/Cora/nn/prop_2.vnnlib'];
[nn2,opts2,x2,r2,A2,b2,safeSet2] = aux_readNetworkAndOptions( ...
    model1Path,prop2Path);

% =========================================================================
% Test each option one at a time against the baseline
% =========================================================================

% --- falsification_method ------------------------------------------------
falsifMethods = {'center','fgsm','zonotack'};
for i = 1:length(falsifMethods)
    opts = opts2;
    opts.nn.falsification_method = falsifMethods{i};
    [verifRes,x_,y_] = nn2.verify(x2,r2,A2,b2,safeSet2,opts,timeout,verbose);
    assertLoop(~strcmp(verifRes.str,'VERIFIED'),i);
    if strcmp(verifRes.str,'COUNTEREXAMPLE')
        aux_validateCounterexample(nn2,x2,r2,A2,b2,safeSet2,x_,y_);
    end
end

% --- refinement_method ---------------------------------------------------
refMethods = {'naive','zonotack'};
for i = 1:length(refMethods)
    opts = opts1;
    opts.nn.refinement_method = refMethods{i};
    [verifRes,x_,y_] = nn1.verify(x1,r1,A1,b1,safeSet1,opts,timeout,verbose);
    assertLoop(~strcmp(verifRes.str,'COUNTEREXAMPLE'),i);
    assertLoop(isempty(x_) & isempty(y_),i);
end

% --- poly_method ---------------------------------------------------------
polyMethods = {'bounds','singh'};
for i = 1:length(polyMethods)
    opts = opts1;
    opts.nn.poly_method = polyMethods{i};
    [verifRes,x_,y_] = nn1.verify(x1,r1,A1,b1,safeSet1,opts,timeout,verbose);
    assertLoop(~strcmp(verifRes.str,'COUNTEREXAMPLE'),i);
    assertLoop(isempty(x_) & isempty(y_),i);
end

% --- num_splits ----------------------------------------------------------
numSplitsVals = [1, 2, 3];
for i = 1:length(numSplitsVals)
    opts = opts1;
    opts.nn.num_splits = numSplitsVals(i);
    [verifRes,x_,y_] = nn1.verify(x1,r1,A1,b1,safeSet1,opts,timeout,verbose);
    assertLoop(~strcmp(verifRes.str,'COUNTEREXAMPLE'),i);
    assertLoop(isempty(x_) & isempty(y_),i);
end

% --- num_dimensions ------------------------------------------------------
numDimsVals = [0, 1, 2];
for i = 1:length(numDimsVals)
    opts = opts1;
    opts.nn.num_dimensions = numDimsVals(i);
    [verifRes,x_,y_] = nn1.verify(x1,r1,A1,b1,safeSet1,opts,timeout,verbose);
    assertLoop(~strcmp(verifRes.str,'COUNTEREXAMPLE'),i);
    assertLoop(isempty(x_) & isempty(y_),i);
end

% --- num_neuron_splits ---------------------------------------------------
neuronSplitVals = [0, 1];
for i = 1:length(neuronSplitVals)
    opts = opts1;
    opts.nn.num_neuron_splits = neuronSplitVals(i);
    [verifRes,x_,y_] = nn1.verify(x1,r1,A1,b1,safeSet1,opts,timeout,verbose);
    assertLoop(~strcmp(verifRes.str,'COUNTEREXAMPLE'),i);
    assertLoop(isempty(x_) & isempty(y_),i);
end

% --- num_relu_constraints ------------------------------------------------
reluConstVals = [0, 3, inf];
for i = 1:length(reluConstVals)
    opts = opts1;
    opts.nn.num_relu_constraints = reluConstVals(i);
    [verifRes,x_,y_] = nn1.verify(x1,r1,A1,b1,safeSet1,opts,timeout,verbose);
    assertLoop(~strcmp(verifRes.str,'COUNTEREXAMPLE'),i);
    assertLoop(isempty(x_) & isempty(y_),i);
end

% --- approx_error_order --------------------------------------------------
approxOrders = {'length','sensitivity*length'};
for i = 1:length(approxOrders)
    opts = opts1;
    opts.nn.approx_error_order = approxOrders{i};
    [verifRes,x_,y_] = nn1.verify(x1,r1,A1,b1,safeSet1,opts,timeout,verbose);
    assertLoop(~strcmp(verifRes.str,'COUNTEREXAMPLE'),i);
    assertLoop(isempty(x_) & isempty(y_),i);
end

% --- train.num_init_gens -------------------------------------------------
initGensVals = [5, inf];
for i = 1:length(initGensVals)
    opts = opts1;
    opts.nn.train.num_init_gens = initGensVals(i);
    [verifRes,x_,y_] = nn1.verify(x1,r1,A1,b1,safeSet1,opts,timeout,verbose);
    assertLoop(~strcmp(verifRes.str,'COUNTEREXAMPLE'),i);
    assertLoop(isempty(x_) & isempty(y_),i);
end

% --- train.num_approx_err ------------------------------------------------
approxErrVals = [25, inf];
for i = 1:length(approxErrVals)
    opts = opts1;
    opts.nn.train.num_approx_err = approxErrVals(i);
    [verifRes,x_,y_] = nn1.verify(x1,r1,A1,b1,safeSet1,opts,timeout,verbose);
    assertLoop(~strcmp(verifRes.str,'COUNTEREXAMPLE'),i);
    assertLoop(isempty(x_) & isempty(y_),i);
end

% --- interval_center -----------------------------------------------------
for i = 1:2
    vals = {true, false};
    opts = opts1;
    opts.nn.interval_center = vals{i};
    [verifRes,x_,y_] = nn1.verify(x1,r1,A1,b1,safeSet1,opts,timeout,verbose);
    assertLoop(~strcmp(verifRes.str,'COUNTEREXAMPLE'),i);
    assertLoop(isempty(x_) & isempty(y_),i);
end

% --- verify_dequeue_type -------------------------------------------------
deqTypes = {'front','half-half'};
for i = 1:length(deqTypes)
    opts = opts1;
    opts.nn.verify_dequeue_type = deqTypes{i};
    [verifRes,x_,y_] = nn1.verify(x1,r1,A1,b1,safeSet1,opts,timeout,verbose);
    assertLoop(~strcmp(verifRes.str,'COUNTEREXAMPLE'),i);
    assertLoop(isempty(x_) & isempty(y_),i);
end

% --- verify_enqueue_type -------------------------------------------------
enqTypes = {'append','prepend'};
for i = 1:length(enqTypes)
    opts = opts1;
    opts.nn.verify_enqueue_type = enqTypes{i};
    [verifRes,x_,y_] = nn1.verify(x1,r1,A1,b1,safeSet1,opts,timeout,verbose);
    assertLoop(~strcmp(verifRes.str,'COUNTEREXAMPLE'),i);
    assertLoop(isempty(x_) & isempty(y_),i);
end

% --- max_verif_iter ------------------------------------------------------
maxIterVals = [1, 5];
for i = 1:length(maxIterVals)
    opts = opts1;
    opts.nn.max_verif_iter = maxIterVals(i);
    [verifRes,x_,y_] = nn1.verify(x1,r1,A1,b1,safeSet1,opts,timeout,verbose);
    assertLoop(~strcmp(verifRes.str,'COUNTEREXAMPLE'),i);
    assertLoop(isempty(x_) & isempty(y_),i);
end

% --- input_generator_heuristic -------------------------------------------
inGenHeurs = {'most-sensitive-input-radius','zono-norm-gradient'};
for i = 1:length(inGenHeurs)
    opts = opts1;
    opts.nn.input_generator_heuristic = inGenHeurs{i};
    [verifRes,x_,y_] = nn1.verify(x1,r1,A1,b1,safeSet1,opts,timeout,verbose);
    assertLoop(~strcmp(verifRes.str,'COUNTEREXAMPLE'),i);
    assertLoop(isempty(x_) & isempty(y_),i);
end

% --- input_split_heuristic -----------------------------------------------
inSplitHeurs = {'most-sensitive-input-radius','zono-norm-gradient'};
for i = 1:length(inSplitHeurs)
    opts = opts1;
    opts.nn.input_split_heuristic = inSplitHeurs{i};
    [verifRes,x_,y_] = nn1.verify(x1,r1,A1,b1,safeSet1,opts,timeout,verbose);
    assertLoop(~strcmp(verifRes.str,'COUNTEREXAMPLE'),i);
    assertLoop(isempty(x_) & isempty(y_),i);
end

% --- neuron_split_heuristic (requires num_neuron_splits >= 1) ------------
neuronHeurs = {'least-unstable','least-unstable-gradient', ...
    'most-sensitive-approx-error','most-sensitive-input-radius', ...
    'zono-norm-gradient'};
for i = 1:length(neuronHeurs)
    opts = opts1;
    opts.nn.num_neuron_splits = 1;
    opts.nn.neuron_split_heuristic = neuronHeurs{i};
    [verifRes,x_,y_] = nn1.verify(x1,r1,A1,b1,safeSet1,opts,timeout,verbose);
    assertLoop(~strcmp(verifRes.str,'COUNTEREXAMPLE'),i);
    assertLoop(isempty(x_) & isempty(y_),i);
end

% --- neuron_split_position (requires num_neuron_splits >= 1) ------------
neuronSplitPos = {'zero','middle'};
for i = 1:length(neuronSplitPos)
    opts = opts1;
    opts.nn.num_neuron_splits = 1;
    opts.nn.neuron_split_position = neuronSplitPos{i};
    [verifRes,x_,y_] = nn1.verify(x1,r1,A1,b1,safeSet1,opts,timeout,verbose);
    assertLoop(~strcmp(verifRes.str,'COUNTEREXAMPLE'),i);
    assertLoop(isempty(x_) & isempty(y_),i);
end

% --- relu_constraint_heuristic (requires num_relu_constraints >= 1) ------
reluHeurs = {'least-unstable','least-unstable-gradient', ...
    'most-sensitive-approx-error','most-sensitive-input-radius', ...
    'zono-norm-gradient'};
for i = 1:length(reluHeurs)
    opts = opts1;
    opts.nn.num_relu_constraints = 3;
    opts.nn.relu_constraint_heuristic = reluHeurs{i};
    [verifRes,x_,y_] = nn1.verify(x1,r1,A1,b1,safeSet1,opts,timeout,verbose);
    assertLoop(~strcmp(verifRes.str,'COUNTEREXAMPLE'),i);
    assertLoop(isempty(x_) & isempty(y_),i);
end

end


% Auxiliary functions -----------------------------------------------------

function [nn,options,x,r,A,b,safeSet] = ...
        aux_readNetworkAndOptions(modelPath,vnnlibPath)
  % Create evaluation options.
  options.nn = struct(...
      'use_approx_error',true,...
      'poly_method','bounds',...
      'train',struct(...
          'backprop',false,...
          'mini_batch_size',2^8 ...
      ) ...
  );
  % Set default training parameters.
  options = nnHelper.validateNNoptions(options,true);

  % Use zonotack baseline.
  options.nn.falsification_method = 'zonotack';
  options.nn.refinement_method = 'zonotack';
  options.nn.num_splits = 2;
  options.nn.num_dimensions = 1;
  options.nn.num_neuron_splits = 0;

  % Read the neural network.
  nn = neuralNetwork.readONNXNetwork(modelPath,false,'BSSC');

  % Read the input set and specification.
  [X0,specs] = vnnlib2cora(vnnlibPath);

  % Extract input set.
  x = 1/2*(X0{1}.sup + X0{1}.inf);
  r = 1/2*(X0{1}.sup - X0{1}.inf);

  % Extract specification.
  if isa(specs.set,'halfspace')
      A = specs.set.c';
      b = specs.set.d;
  else
      A = specs.set.A;
      b = specs.set.b;
  end
  safeSet = strcmp(specs.type,'safeSet');
end

function aux_validateCounterexample(nn,x,r,A,b,safeSet,x_,y_)
    assert(~isempty(x_) & ~isempty(y_));
    % Check input bounds.
    tol = 1e-6;
    assert(all(x_ >= x - r - tol,'all') & all(x_ <= x + r + tol,'all'), ...
        'Counterexample x_ out of input bounds.');
    % Check output matches.
    yi = nn.evaluate(x_);
    assert(all(abs(y_ - yi) <= 1e-7,'all'), ...
        'Counterexample y_ does not match nn.evaluate(x_).');
    % Check spec is violated.
    if safeSet
        assert(any(A*yi >= b,1), ...
            'Counterexample does not violate safe-set specification.');
    else
        assert(all(A*yi <= b,1), ...
            'Counterexample does not violate unsafe-set specification.');
    end
end

% ------------------------------ END OF CODE ------------------------------
