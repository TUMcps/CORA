function res = testnn_verify()
% testnn_verify - test neuralNetwork/verify function
%
% Syntax:
%    res = testnn_verify()
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

modelPath = [CORAROOT '/models/Cora/nn/ACASXU_run2a_1_2_batch_2000.onnx'];
timeout = 100;

% First test case: prop_1.vnnlib
prop1Filename = [CORAROOT '/models/Cora/nn/prop_1.vnnlib'];
[nn,options,x,r,A,b,safeSet] = ...
    aux_readNetworkAndOptions(modelPath,prop1Filename);
% Do verification.
[verifRes,x_,y_] = nn.verify(x,r,A,b,safeSet,options,timeout,false);
assert(strcmp(verifRes,'VERIFIED') & isempty(x_) & isempty(y_));

% Second test case: prop_2.vnnlib
prop2Filename = [CORAROOT '/models/Cora/nn/prop_2.vnnlib'];
[nn,options,x,r,A,b,safeSet] = ...
    aux_readNetworkAndOptions(modelPath,prop2Filename);
% Do verification.
[verifRes,x_,y_] = nn.verify(x,r,A,b,safeSet,options,timeout,false);
assert(strcmp(verifRes,'COUNTEREXAMPLE') & ~isempty(x_) & ~isempty(y_));
assert(aux_checkCounterexample(nn,A,b,safeSet,x_,y_));

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
          'mini_batch_size',512 ...
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
      b = -specs.set.d;
  else
      A = specs.set.A;
      b = -specs.set.b;
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
    violates = any(A*yi + b >= 0,1);
else
    violates = all(A*yi + b <= 0,1);
end
assert(res & violates);
end

% ------------------------------ END OF CODE ------------------------------
