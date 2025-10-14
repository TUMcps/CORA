function res = testnn_neuralNetwork_verify()
% testnn_neuralNetwork_verify - test neuralNetwork/verify function
%
% Syntax:
%    res = testnn_neuralNetwork_verify()
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

% Authors:       Lukas Koller
% Written:       03-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% We use the specs from the acasxu benchmark: prop_1, prop_2, prop_3, and
% prop_5.

res = true;

% Toggle verbose verification output.
verbose = true;

% Specify the model path.
model1Path = [CORAROOT '/models/Cora/nn/ACASXU_run2a_1_2_batch_2000.onnx'];
model2Path = [CORAROOT '/models/Cora/nn/ACASXU_run2a_5_3_batch_2000.onnx'];
prop1Filename = [CORAROOT '/models/Cora/nn/prop_1.vnnlib'];
prop2Filename = [CORAROOT '/models/Cora/nn/prop_2.vnnlib'];

% Set a timeout of 2s.
timeout = 10;

% First test case: prop_1.vnnlib ------------------------------------------
[nn,options,x,r,A,b,safeSet] = ...
    aux_readNetworkAndOptions(model1Path,prop1Filename);
% Test 'naive'-splitting and 'fgsm'-falsification.
options.nn.falsification_method = 'fgsm';
options.nn.refinement_method = 'naive';
% Do verification.
[verifRes,x_,y_] = nn.verify(x,r,A,b,safeSet,options,timeout,verbose);
assert(~strcmp(verifRes.str,'COUNTEREXAMPLE') & isempty(x_) & isempty(y_));
% assert(strcmp(verifRes.str,'VERIFIED') & isempty(x_) & isempty(y_));

% Test 'zonotack' implementation with restricted number of generators.
options.nn.falsification_method = 'zonotack';
options.nn.refinement_method = 'zonotack';
% Specify parameters.
options.nn.num_splits = 2; 
options.nn.num_dimensions = 2;
options.nn.num_neuron_splits = 1;
% Restrict the number of input generators.
options.nn.train.num_init_gens = 5;
% Restrict the number of approximation error generators per layer.
options.nn.train.num_approx_err = 50;
options.nn.approx_error_order = 'sensitivity*length';
% Add relu tightening constraints.
options.nn.num_relu_tighten_constraints = inf;
% options.nn.neuron_xor_input_splits = false;
% Do verification.
[verifRes,x_,y_] = nn.verify(x,r,A,b,safeSet,options,timeout,verbose);
assert(~strcmp(verifRes.str,'COUNTEREXAMPLE') & isempty(x_) & isempty(y_));
% assert(strcmp(verifRes.str,'VERIFIED') & isempty(x_) & isempty(y_));

% Second test case: prop_2.vnnlib -----------------------------------------
[nn,options,x,r,A,b,safeSet] = ...
    aux_readNetworkAndOptions(model1Path,prop2Filename);
% Test 'naive'-splitting and 'fgsm'-falsification.
options.nn.falsification_method = 'fgsm';
options.nn.refinement_method = 'naive';
% Do verification.
[verifRes,x_,y_] = nn.verify(x,r,A,b,safeSet,options,timeout,verbose);
assert(~strcmp(verifRes.str,'VERIFIED'));
if strcmp(verifRes.str,'COUNTEREXAMPLE')
    assert(~isempty(x_) & ~isempty(y_) & ...
        aux_checkCounterexample(nn,A,b,safeSet,x_,y_));
end

% Test 'zonotack' implementation with restricted number of generators.
options.nn.falsification_method = 'zonotack';
options.nn.refinement_method = 'zonotack-layerwise';
% Specify parameters.
options.nn.num_splits = 3; 
options.nn.num_dimensions = 1;
options.nn.num_neuron_splits = 0;
% Restrict the number of input generators.
options.nn.train.num_init_gens = 5;
% Restrict the number of approximation error generators per layer.
options.nn.train.num_approx_err = 25;
options.nn.approx_error_order = 'length';
% Add relu tightening constraints.
options.nn.num_relu_tighten_constraints = 3;
% Do verification.
[verifRes,x_,y_] = nn.verify(x,r,A,b,safeSet,options,timeout,verbose);
assert(~strcmp(verifRes.str,'VERIFIED'));
if strcmp(verifRes.str,'COUNTEREXAMPLE')
    assert(~isempty(x_) & ~isempty(y_) & ...
        aux_checkCounterexample(nn,A,b,safeSet,x_,y_));
end

% Test 'zonotack' implementation with restricted number of generators.
options.nn.falsification_method = 'zonotack';
options.nn.refinement_method = 'zonotack-layerwise';
% Specify parameters.
options.nn.num_splits = 2; 
options.nn.num_dimensions = 1;
options.nn.num_neuron_splits = 0;
% Add relu tightening constraints.
options.nn.num_relu_tighten_constraints = 100;
% Do verification.
[verifRes,x_,y_] = nn.verify(x,r,A,b,safeSet,options,timeout,verbose);
assert(~strcmp(verifRes.str,'VERIFIED'));
if strcmp(verifRes.str,'COUNTEREXAMPLE')
    assert(~isempty(x_) & ~isempty(y_) & ...
        aux_checkCounterexample(nn,A,b,safeSet,x_,y_));
end

% Third test case with other model: prop_2.vnnlib -------------------------
% Specify the model path.
[nn,options,x,r,A,b,safeSet] = ...
    aux_readNetworkAndOptions(model2Path,prop2Filename);
options.nn.falsification_method = 'zonotack';
options.nn.refinement_method = 'zonotack'; 
% Specify parameters.
options.nn.num_splits = 2; 
options.nn.num_dimensions = 0;
options.nn.num_neuron_splits = 1;
options.nn.add_orth_neuron_splits = true;
% Add relu tightening constraints.
options.nn.num_relu_constraints = 15;
% Do verification.
[verifRes,x_,y_] = nn.verify(x,r,A,b,safeSet,options,timeout,verbose);
% Finding a counterexample is hard.
assert(~strcmp(verifRes.str,'VERIFIED'));
if strcmp(verifRes.str,'COUNTEREXAMPLE')
    assert(~isempty(x_) & ~isempty(y_) & ...
        aux_checkCounterexample(nn,A,b,safeSet,x_,y_));
end

% Fourth test case with other model: prop_2.vnnlib -------------------------
% Specify the model path.
[nn,options,x,r,A,b,safeSet] = ...
    aux_readNetworkAndOptions(model2Path,prop2Filename);
options.nn.falsification_method = 'zonotack';
options.nn.refinement_method = 'zonotack'; 
% Specify parameters.
options.nn.num_splits = 2; 
options.nn.num_dimensions = 1;
options.nn.num_neuron_splits = 1;
options.nn.input_xor_neuron_splitting = true;
% Add relu tightening constraints.
options.nn.num_relu_constraints = inf;
% Do verification.
[verifRes,x_,y_] = nn.verify(x,r,A,b,safeSet,options,timeout,verbose);
% Finding a counterexample is hard.
assert(~strcmp(verifRes.str,'VERIFIED'));
if strcmp(verifRes.str,'COUNTEREXAMPLE')
    assert(~isempty(x_) & ~isempty(y_) & ...
        aux_checkCounterexample(nn,A,b,safeSet,x_,y_));
end

end


% Auxiliary functions -----------------------------------------------------

function [nn,options,x,r,A,b,safeSet] = ...
        aux_readNetworkAndOptions(modelPath,vnnlibPath)
  % Create evaluation options.
  options.nn = struct(...
      'use_approx_error',true,...
      'poly_method','bounds',...'bounds','singh'
      'train',struct(...
          'backprop',false,...
          'mini_batch_size',2^8 ...
      ) ...
  );
  % Set default training parameters
  options = nnHelper.validateNNoptions(options,true);
  options.nn.interval_center = false;

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

function res = aux_checkCounterexample(nn,A,b,safeSet,x_,y_)
% Compute output of the neural network.
yi = nn.evaluate(x_);
% Check if output matches.
res = all(abs(y_ - yi) <= 1e-7,'all');
% Check of output violates the specification.
if safeSet
    violates = any(A*yi >= b,1);
else
    violates = all(A*yi <= b,1);
end
assert(res & violates);
end

% ------------------------------ END OF CODE ------------------------------
