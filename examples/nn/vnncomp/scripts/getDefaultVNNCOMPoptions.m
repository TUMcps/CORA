function options = getDefaultVNNCOMPoptions(benchName)
% getDefaultVNNCOMPoptions - return the default verification options for a benchmark.
%
% Syntax:
%    options = getDefaultVNNCOMPoptions(benchName)
%
% Inputs:
%    benchName - name of the benchmark (e.g. 'acasxu_2023', 'safenlp_2024')
%
% Outputs:
%    options - struct with options.nn set to benchmark defaults
%
% See also: prepare_instance, hyperparam_tuning_vnncomp

% Authors:       Benedikt Kellner, Lukas Koller
% Written:       04-May-2026
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

options.nn = struct(...
    'use_approx_error',true,...
    'poly_method','bounds',... % {'bounds','singh','center'}
    'train',struct(...
        'backprop',false,...
        'mini_batch_size',2^10 ...
    ) ...
);
options = nnHelper.validateNNoptions(options,true);
options.nn.interval_center = false;
options.nn.batch_norm_moving_stats = true;
options.nn.falsification_method = 'zonotack'; % {'center','fgsm','zonotack'}
options.nn.refinement_method = 'zonotack';    % {'naive','zonotack','zonotack-layerwise'}
options.nn.train.num_init_gens = inf;
options.nn.approx_error_order = 'sensitivity*length';
options.nn.conzonotope_bounding_method = 'dual-iter';
options.nn.num_splits = 2;
options.nn.num_dimensions = 1;
options.nn.num_neuron_splits = 0;
options.nn.num_relu_constraints = 0;
options.nn.input_xor_neuron_splitting = true;
options.nn.polytope_bound_approx_max_iter = 3;
options.nn.refinement_min_iter = 4;
options.nn.refinement_max_iter = 8;
options.nn.input_generator_heuristic = 'zono-norm-gradient';
options.nn.input_split_heuristic = 'zono-norm-gradient';
options.nn.neuron_split_heuristic = 'zono-norm-gradient';
options.nn.relu_constraint_heuristic = 'zono-norm-gradient';

% strip year suffix, e.g. 'acasxu_2023' -> 'acasxu'
tok = regexp(benchName,'^(.+?)[_-](20\d{2})$','tokens');
if ~isempty(tok)
    benchName_ = tok{1}{1};
else
    benchName_ = benchName;
end

% benchmark-specific option overrides (VNN-COMP'25)
if strcmp(benchName_,'cgan')
    options.nn.use_dlconv = true;
    options.nn.num_splits = 2;
    options.nn.num_dimensions = 1;
    options.nn.num_neuron_splits = 1;
    options.nn.train.mini_batch_size = 2^2;
    options.nn.max_verif_iter = 100;
    options.nn.neuron_split_heuristic = 'least-unstable';
    options.nn.exact_conzonotope_bounds = true;
    options.nn.verify_cascade_unsafe_set_constraints = false;
    options.nn.num_relu_constraints = 0;

elseif strcmp(benchName_,'cifar100') % large image classification
    options.nn.train.backprop = true;  % required for interval gradient (input generator heuristic)
    options.nn.train.mini_batch_size = 8;
    options.nn.input_split_heuristic = 'zono-norm-gradient';
    options.nn.neuron_split_heuristic = 'least-unstable';
    options.nn.refinement_method = 'zonotack';
    options.nn.falsification_method = 'fgsm';
    options.nn.num_split_dimensions = 0;
    options.nn.num_splits = 2;
    options.nn.num_neuron_splits = 1;
    options.nn.input_xor_neuron_splitting = false;
    options.nn.num_relu_constraints = 100;
    options.nn.train.num_init_gens = 250;
    options.nn.train.num_approx_err = 25;
    options.nn.exact_conzonotope_bounds = false;
    options.nn.max_verif_iter = 1000;
    options.nn.verify_dequeue_type = 'half-half';
    options.nn.verify_enqueue_type = 'append';

elseif strcmp(benchName_,'collins_rul_cnn') % VNN-COMP'24
    options.nn.interval_center = true;
    options.nn.train.num_init_gens = inf;
    options.nn.train.num_approx_err = 100;

elseif strcmp(benchName_,'cora')
    options.nn.num_relu_constraints = inf;
    options.nn.train.mini_batch_size = 2^5;
    options.nn.num_splits = 2;
    options.nn.num_dimensions = 1;
    options.nn.num_neuron_splits = 1;
    options.nn.batch_union_conzonotope_bounds = false;

elseif strcmp(benchName_,'metaroom')
    options.nn.train.num_init_gens = 500;
    options.nn.train.num_approx_err = 100;
    options.nn.train.mini_batch_size = 2^2;
    options.nn.num_splits = 2;
    options.nn.num_dimensions = 1;
    options.nn.num_neuron_splits = 0;
    options.nn.num_relu_constraints = 0;
    options.nn.batch_union_conzonotope_bounds = false;
    options.nn.max_verif_iter = 10;

elseif strcmp(benchName_,'mnist_fc') % VNN-COMP'22
    options.nn.interval_center = true;
    options.nn.train.num_init_gens = 500;
    options.nn.train.num_approx_err = 100;
    options.nn.batch_union_conzonotope_bounds = false;

elseif strcmp(benchName_,'oval21')
    options.nn.interval_center = true;
    options.nn.train.num_init_gens = 100;
    options.nn.train.num_approx_err = 10;
    options.nn.batch_union_conzonotope_bounds = false;

elseif ismember(benchName_,{'monotonic_acasxu','isomorphic_acasxu'})
    % VNN-COMP'26 multi-network ACAS Xu benchmarks.
    % Small composite networks (joint f+g, 5-in each); use input-radius
    % heuristics to avoid backprop through the composite layer.
    options.nn.num_splits = 2;
    options.nn.num_dimensions = 1;
    options.nn.num_neuron_splits = 0;
    options.nn.train.mini_batch_size = 2^10;

elseif strcmp(benchName_,'safenlp')
    options.nn.num_splits = 2;
    options.nn.num_dimensions = 1;
    options.nn.num_neuron_splits = 1;
    options.nn.num_relu_constraints = 100;

elseif strcmp(benchName_,'tinyimagenet') % VNN-COMP'24
    options.nn.interval_center = true;
    options.nn.train.num_init_gens = 500;
    options.nn.train.num_approx_err = 0;
    options.nn.num_relu_constraints = 0;
    options.nn.num_splits = 2;
    options.nn.num_dimensions = 1;
    options.nn.num_neuron_splits = 1;
    options.nn.batch_union_conzonotope_bounds = false;
    options.nn.train.mini_batch_size = 2^5;

end

end

% ------------------------------ END OF CODE ------------------------------
