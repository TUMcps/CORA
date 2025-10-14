function res = prepare_instance(benchName,modelPath,vnnlibPath)
% prepare_instance - called before the verification timer is started. Load
% the neural network and the vnnlib specification; they are stored as a
% .mat file, which is opened by run_instance.m
%
% Syntax:
%    res = prepare_instance(benchName,modelPath,vnnlibPath)
%
% Inputs:
%    benchName - name of the benchmark
%    modelPath - path to the .onnx-file
%    vnnlibPath - path to the .vnnlib-file
%
% Outputs:
%    res - result code
%
% References:
%    [1] VNN-COMP'24
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Lukas Koller
% Written:       11-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

  fprintf('prepare_instance(%s,%s,%s)...\n',benchName,modelPath,vnnlibPath);
  try
      fprintf('--- Loading network...');
      % Load neural network.
      [nn,options,permuteDims] = aux_readNetworkAndOptions( ...
          benchName,modelPath,vnnlibPath,false);

      fprintf(' done\n');

      fprintf('--- GPU available: %d\n',options.nn.train.use_gpu);
   
      fprintf('--- Loading specification...');
      % Load specification.
      [X0,specs] = vnnlib2cora(vnnlibPath);
      fprintf(' done\n');

      fprintf('--- Storing MATLAB file...');
      % Create filename.
      instanceFilename = getInstanceFilename(benchName,modelPath,vnnlibPath);
      % Store network, options, and specification.
      save(instanceFilename,'nn','options','permuteDims','X0','specs');
      fprintf(' done\n');

      % Print the options.
      aux_printOptions(options);
  catch e
      % Print the error message. 
      fprintf(newline);
      fprintf(e.message);
      fprintf(newline);
      % Some error
      res = 1;
      return;
  end
  res = 0;
end


% Auxiliary functions -----------------------------------------------------

function [nn,options,permuteDims] = aux_readNetworkAndOptions( ...
  benchName,modelPath,vnnlibPath,verbose)

  % Create evaluation options.
  options.nn = struct(...
      'use_approx_error',true,...
      'poly_method','bounds',... {'bounds','singh','center'}
      'train',struct(...
          ... 'use_gpu',false,...
          'backprop',false,...
          'mini_batch_size',2^10 ...
      ) ...
  );
  % Set default training parameters
  options = nnHelper.validateNNoptions(options,true);
  % Disable the interval-center by default.
  options.nn.interval_center = false;
  % Use the moving statistics for the batch normalization.
  options.nn.batch_norm_moving_stats = true;

  % Specify falsification method: {'center','fgsm','zonotack'}.
  options.nn.falsification_method = 'zonotack';
  % Specify input set refinement method: {'naive','zonotack','zonotack-layerwise'}.
  options.nn.refinement_method = 'zonotack';
  % Set number of input generators.
  options.nn.train.num_init_gens = inf;
  % Set number of approximation error generators per layer.
  options.nn.approx_error_order = 'length'; % 'sensitivity*'
  % Compute the exact bounds of the constraint zonotope.
  options.nn.exact_conzonotope_bounds = false;
  % Specify number of splits, dimensions, and neuron-splits.
  options.nn.num_splits = 2; 
  options.nn.num_dimensions = 1;
  options.nn.num_neuron_splits = 0;
  % Add relu tightening constraints.
  options.nn.num_relu_constraints = 0;
  % Add orthogonal split constraints.
  options.nn.add_orth_neuron_splits = false;
  % Either neuron or input splitting.
  options.nn.input_xor_neuron_splitting = true;
  % Specify the number of iterations.
  options.nn.polytope_bound_approx_max_iter = 3;
  options.nn.refinement_max_iter = 3;
  % Specify the queue-style.
  % options.nn.verify_dequeue_type = 'half-half';
  % options.nn.verify_enqueue_type = 'prepend';

  % Default: do not permute the input dimensions. 
  permuteDims = false;

  % Obtain the model name.
  [~,modelName,~] = getInstanceFilename(benchName,modelPath,vnnlibPath);

  % VNN-COMP'24 Benchmarks ------------------------------------------------
  if strcmp(benchName,'test') ...
    || strcmp(modelName{1},'test_nano') % is called after each benchmark
      nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'','', ...
          'dlnetwork',false);
      % Use the default parameters.
  elseif strcmp(benchName,'acasxu_2023')
    % acasxu ----------------------------------------------------------
    nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'BSSC');
    % Specify number of splits, dimensions, and neuron-splits.
    options.nn.num_splits = 2; 
    options.nn.num_dimensions = 1;
    % options.nn.num_neuron_splits = 1;
    % Add relu tightening constraints.
    % options.nn.num_relu_constraints = inf;
  elseif strcmp(benchName,'cctsdb_yolo_2023')
      throw(CORAerror('CORA:notSupported',...
          sprintf("Benchmark '%s' not supported!",benchName)));
  elseif strcmp(benchName,'cgan_2023')
      % c_gan -----------------------------------------------------------
      % nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'BC');
      % --- TODO: implement convTranspose
      throw(CORAerror('CORA:notSupported',...
          sprintf("Benchmark '%s' not supported!",benchName)));
  elseif strcmp(benchName,'cifar100')
      % vnncomp2024_cifar100_benchmark ----------------------------------
      nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'BCSS', ...
          '','dagnetwork',true);
      % Bring input into the correct shape.
      permuteDims = true;
      % Requires less memory.
      % options.nn.falsification_method = 'fgsm';
      % Use interval-center.
      options.nn.interval_center = true;
      options.nn.train.num_init_gens = 1000;
      options.nn.train.num_approx_err = 100;
      % Add relu tightening constraints.
      % options.nn.num_relu_constraints = 100;
      % Specify number of splits, dimensions, and neuron-splits.
      options.nn.num_splits = 2; 
      options.nn.num_dimensions = 1;
      options.nn.num_neuron_splits = 1;
      % Save memory (reduce batch size & do not batch union constraints).
      options.nn.train.mini_batch_size = 2^2;
      options.nn.batch_union_conzonotope_bounds = false;
      % To find all counterexamples we restrict the number of verification
      % iterations.
      options.nn.max_verif_iter = 10;
  elseif strcmp(benchName,'collins_aerospace_benchmark')
      throw(CORAerror('CORA:notSupported',...
          sprintf("Benchmark '%s' not supported!",benchName)));
  elseif strcmp(benchName,'collins_rul_cnn_2023') ...
      || strcmp(benchName,'collins_rul_cnn')
      % collins_rul_cnn -------------------------------------------------
      nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'BCSS');
      % Bring input into the correct shape.
      permuteDims = true;
      % Use interval-center.
      options.nn.interval_center = true;
      options.nn.train.num_init_gens = inf;
      options.nn.train.num_approx_err = 100;
  elseif strcmp(benchName,'collins_yolo_robustness_2023')
      throw(CORAerror('CORA:notSupported',...
          sprintf("Benchmark '%s' not supported!",benchName)));
  elseif strcmp(benchName,'cora')
      nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'BC');
      % Requires less memory.
      % options.nn.falsification_method = 'fgsm';
      % Use the default parameters.
      % options.nn.interval_center = true;
      % options.nn.train.num_init_gens = 500; % inf;
      % options.nn.train.num_approx_err = 100;
      % Add relu tightening constraints.
      options.nn.num_relu_constraints = inf;
      % Reduce batch size.
      options.nn.train.mini_batch_size = 2^5;
      % Specify number of splits, dimensions, and neuron-splits.
      options.nn.num_splits = 2; 
      options.nn.num_dimensions = 1;
      options.nn.num_neuron_splits = 1;
      % Save memory (do not batch union constraints).
      options.nn.batch_union_conzonotope_bounds = false;
  elseif strcmp(benchName,'dist_shift_2023')
      % dist_shift ------------------------------------------------------
      nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'BC');
      % Use default values.
  elseif strcmp(benchName,'linearizenn')
      % LinearizeNN -----------------------------------------------------
      nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'BC', ...
          '','dagnetwork',true);
      % Specify number of splits, dimensions, and neuron-splits.
      options.nn.num_splits = 2; 
      options.nn.num_dimensions = 1;
      options.nn.num_neuron_splits = 0;
  elseif strcmp(benchName,'lsnc')
      throw(CORAerror('CORA:notSupported',...
          sprintf("Benchmark '%s' not supported!",benchName)));
  elseif strcmp(benchName,'metaroom_2023')
      % metaroom --------------------------------------------------------
      nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'BCSS');
      % Bring input into the correct shape.
      permuteDims = true;
      % Use interval-center.
      % options.nn.interval_center = true;
      options.nn.train.num_init_gens = 500;
      options.nn.train.num_approx_err = 100;
      % Reduce the batch size.
      options.nn.train.mini_batch_size = 2^2;
      % Specify number of splits, dimensions, and neuron-splits.
      options.nn.num_splits = 2; 
      options.nn.num_dimensions = 1;
      options.nn.num_neuron_splits = 0;
      % Add relu tightening constraints.
      options.nn.num_relu_constraints = 0; % inf;
      % Save memory (do not batch union constraints).
      options.nn.batch_union_conzonotope_bounds = false;
      % To find all counterexamples we restrict the number of verification
      % iterations.
      options.nn.max_verif_iter = 10;
  elseif strcmp(benchName,'ml4acopf_2023')
      throw(CORAerror('CORA:notSupported',...
          sprintf("Benchmark '%s' not supported!",benchName)));
  elseif strcmp(benchName,'ml4acopf_2024')
      throw(CORAerror('CORA:notSupported',...
          sprintf("Benchmark '%s' not supported!",benchName)));
  elseif strcmp(benchName,'nn4sys_2023')
      % nn4sys ----------------------------------------------------------
      if ~strcmp(modelName{1},'lindex') && ...
              ~strcmp(modelName{1},'lindex_deep')
          % Skip this instance.
          throw(CORAerror('CORA:notSupported',...
              sprintf("Model '%s' of benchmark '%s' is not " + ...
              "supported!",modelPath,benchName)));
      end
      nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'BC','BC');
      % Use the default parameters.
  elseif strcmp(benchName,'safenlp')
      % safeNLP ---------------------------------------------------------
      nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'BC');
      % Increase batch size.
      % options.nn.train.mini_batch_size = 2^9;
      % Specify number of splits, dimensions, and neuron-splits.
      options.nn.num_splits = 2; 
      options.nn.num_dimensions = 1;
      options.nn.num_neuron_splits = 1;
      % Add relu tightening constraints.
      options.nn.num_relu_constraints = 100; % inf;
  elseif strcmp(benchName,'tinyimagenet')
      % tinyimagenet ----------------------------------------------------
      nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'BCSS', ...
          '','dagnetwork',true);
      % Bring input into the correct shape.
      permuteDims = true;
      % Requires less memory.
      % options.nn.falsification_method = 'fgsm';
      % Use interval-center.
      options.nn.interval_center = true;
      options.nn.train.num_init_gens = 500;
      options.nn.train.num_approx_err = 0;
      % Add relu tightening constraints.
      options.nn.num_relu_constraints = 0;
      % Specify number of splits, dimensions, and neuron-splits.
      options.nn.num_splits = 2; 
      options.nn.num_dimensions = 1;
      options.nn.num_neuron_splits = 1;
      % Save memory (reduce batch size & do not batch union constraints).
      options.nn.batch_union_conzonotope_bounds = false;
      % Reduce the batch size.
      options.nn.train.mini_batch_size = 2^5;
  elseif strcmp(benchName,'tllverifybench_2023')
      % tllverifybench --------------------------------------------------
      nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'BC');
      % Specify number of splits, dimensions, and neuron-splits.
      % options.nn.num_splits = 2; 
      % options.nn.num_dimensions = 1;
      % options.nn.num_neuron_splits = 1;
  elseif strcmp(benchName,'traffic_signs_recognition_2023')
      throw(CORAerror('CORA:notSupported',...
          sprintf("Benchmark '%s' not supported!",benchName)));
  elseif strcmp(benchName,'vggnet16_2023')
      throw(CORAerror('CORA:notSupported',...
          sprintf("Benchmark '%s' not supported!",benchName)));
  elseif strcmp(benchName,'vit_2023')
      throw(CORAerror('CORA:notSupported',...
          sprintf("Benchmark '%s' not supported!",benchName)));
  elseif strcmp(benchName,'yolo_2023')
      throw(CORAerror('CORA:notSupported',...
          sprintf("Benchmark '%s' not supported!",benchName)));
  % VNN-COMP'22 Benchmarks ------------------------------------------------
  elseif strcmp(benchName,'mnist_fc')
      nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'SSC');
      % Use interval-center.
      options.nn.interval_center = true;
      options.nn.train.num_init_gens = 500;
      options.nn.train.num_approx_err = 100;
      % Save memory (do not batch union constraints).
      options.nn.batch_union_conzonotope_bounds = false;
  elseif strcmp(benchName,'oval21')
      nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'BCSS');
      % Bring input into the correct shape.
      permuteDims = true;
      % Use interval-center.
      options.nn.interval_center = true;
      options.nn.train.num_init_gens = 100;
      options.nn.train.num_approx_err = 10;
      % Save memory (do not batch union constraints).
      options.nn.batch_union_conzonotope_bounds = false;
  elseif strcmp(benchName,'reach_prob_density')
      nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'BC');
  elseif strcmp(benchName,'rl_benchmarks')
      nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'BC');

      if ismember({vnnlibPath}, ...
              {'vnnlib/dubinsrejoin_case_safe_10.vnnlib', ...
              'vnnlib/dubinsrejoin_case_safe_13.vnnlib', ...
              'vnnlib/dubinsrejoin_case_safe_15.vnnlib', ...
              'vnnlib/dubinsrejoin_case_safe_16.vnnlib', ...
              'vnnlib/dubinsrejoin_case_safe_17.vnnlib'})
          % Skip these instances (TODO: fix vnnlib2cora).
          throw(CORAerror('CORA:notSupported',...
              sprintf("Specification '%s' of benchmark '%s' is not " + ...
              "supported!",vnnlibPath,benchName)));
      end
  elseif strcmp(benchName,'sri_resnet_a')
      nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'BCSS');
      % Bring input into the correct shape.
      permuteDims = true;
  elseif strcmp(benchName,'sri_resnet_b')
      nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'BCSS');
      % Bring input into the correct shape.
      permuteDims = true;
  else
      throw(CORAerror('CORA:notSupported',...
          sprintf("Unknown benchmark '%s'!",benchName)));
  end

end

% Auxiliary functions -----------------------------------------------------

function aux_printOptions(options)
    % Print parameters.
    table = CORAtableParameters('neuralNetwork/verify options');
    table.printHeader();
    % Zonotope propagation options.
    table.printContentRow('GPU',string(options.nn.train.use_gpu));
    table.printContentRow('Poly. Method',options.nn.poly_method);
    table.printContentRow('Batchsize', ...
        string(options.nn.train.mini_batch_size));
    table.printContentRow('Interval Center', ...
        string(options.nn.interval_center));
    table.printContentRow('Num. init. Generators', ...
        string(options.nn.train.num_init_gens));
    table.printContentRow('Num. approx. Error (per nonl. Layer)', ...
        string(options.nn.train.num_approx_err));
    table.printContentRow('approx. Error Heuristic', ...
        options.nn.approx_error_order);
    % Main algorithm options.
    table.printContentRow('Falsification Method', ...
        options.nn.falsification_method);
    table.printContentRow('Refinement Method', ...
        options.nn.refinement_method);
    % Details algorithm hyperparameters.
    table.printContentRow('Num. of Splits', ...
        string(options.nn.num_splits));
    table.printContentRow('Num. of Dimensions', ...
        string(options.nn.num_dimensions));
    table.printContentRow('Num. of Neuron-Splits', ...
        string(options.nn.num_neuron_splits));
    table.printContentRow('Num. of ReLU-Tightening Constraints', ...
        string(options.nn.num_relu_constraints));
    % Finish table.
    table.printFooter();
end

% ------------------------------ END OF CODE ------------------------------
