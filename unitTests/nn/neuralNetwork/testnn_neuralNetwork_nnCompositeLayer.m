function res = testnn_neuralNetwork_nnCompositeLayer()
% testnn_neuralNetwork_nnCompositeLayer - parse a neural network from onnx
% containing a composite layer and test its evaluation results.
%
% Syntax:
%    res = testnn_neuralNetwork_nnCompositeLayer()
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
% Written:       27-June-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% We use the networks from cersyve benchmark of the VNN-COMP'25 and check
% evaluation results against the deep learning toolbox.

% Specify the number of test samples.
N = 1000;

% Specify the model path.
model1Path = [CORAROOT '/models/Cora/nn/CERSYVE_cart_pole_finetune_con.onnx']; % Contains 'add'
model2Path = [CORAROOT '/models/Cora/nn/CERSYVE_point_mass_finetune_con.onnx'];
model3Path = [CORAROOT '/models/Cora/nn/LINEARIZENN_AllInOne_10_10.onnx']; % Contains 'concat'
prop1Filename = [CORAROOT '/models/Cora/nn/prop_cart_pole.vnnlib'];
prop2Filename = [CORAROOT '/models/Cora/nn/prop_point_mass.vnnlib'];
prop3Filename = [CORAROOT '/models/Cora/nn/prop_10_10.vnnlib'];

% Load the networks and check their outputs.
res = aux_checkSamples(model1Path,prop1Filename,N);
res = aux_checkSamples(model2Path,prop2Filename,N);
res = aux_checkSamples(model3Path,prop3Filename,N);

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
  nn = neuralNetwork.readONNXNetwork(modelPath,false,'BC','', ...
      'dagnetwork',true);

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

function res = aux_checkSamples(modelPath,propPath,N)
    % Load the networks.
    [nn,~,x,r,~,~,~] = aux_readNetworkAndOptions(modelPath,propPath);
    % Read deep learning toolbox network.
    dlt_nn = importNetworkFromONNX(modelPath);
    
    % Construct the input sets.
    X0 = interval(x - r,x + r);

    % Sample points from the input set.
    xs = randPoint(X0,N);

    % Compute results.
    ys = nn.evaluate(xs);
    
    % Convert the inputs.
    xs_ = dlarray(xs','BC');
    % Initialize the deep learning toolbox network.
    dlt_nn = dlt_nn.initialize(xs_);
    % Compute results.
    ys_ = dlt_nn.predict(xs_);
    try
        % Remove auxiliary directory that is automatically created.
        [~,modelName,~] = fileparts(modelPath);
        rmdir(sprintf('+%s',modelName),'s')
        rmdir('+DLT_CustomLayers','s')
    catch e
    end
    
    % Check results.
    assert(all(withinTol(ys,ys_',1e-6),'all'));

    res = true;
end

% ------------------------------ END OF CODE ------------------------------
