function r = evaluate(obj, input, varargin)
% evaluate - compute the output of a neural network for the given input
%
% Syntax:
%    res = evaluate(obj, input)
%    res = evaluate(obj, input, options)
%
% Inputs:
%    obj - object of class neuralNetwork
%    input - input represented as a numeric or set
%    options - options for neural network evaluation 
%       (stored in options.nn)
%       .bound_approx: bool whether bounds should be overapproximated,
%           or "sample" (not safe!)
%       .reuse_bounds: wheter bounds should be reused
%       .poly_method: how approximation polynomial is found in nonlinear
%           layers, e.g. 'regression', 'singh', 'taylor' ...
%       .num_generators: max number of generators for order reduction
%       .add_approx_error_to_GI: whether 'd' should be added to GI
%       .plot_multi_layer_approx_info: plotting for nnApproximationLayer
%       .max_bounds: max order used in refinement
%       .do_pre_order_reduction: wheter to do a order reduction before
%           evaluating the pZ using the polynomial
%       .max_gens_post: max num_generator before post order reduction
%       .remove_GI: whether to restructue s.t. there remain no ind. gens
%       .force_approx_lin_at: (l, u) distance at which to use 'lin'
%           instead of respective order using 'adaptive'
%       .sort_exponents: whether exponents should be sorted
%       .propagate_bounds: whether bounds should be propagated to next
%           activation layer using interval arithmetic
%       .maxpool_type: for set-based prediction, 'project' or 'regression'
%       .order_reduction_sensitivity: whether sensitivity should be used
%           during order reduction
%       .G: graph object used for graph neural networks
%    idxLayer - indices of layers that should be evaluated
%
% Outputs:
%    res - output of the neural network
%
% References:
%    [1] Kochdumper, N., et al. (2023). Open-and closed-loop neural network
%        verification using polynomial zonotopes. NASA Formal Methods.
%    [2] Ladner, T., et al. (2023). Automatic abstraction refinement in
%        neural network verification using sensitivity analysis. HSCC '23:
%        Proceedings of the 26th International Conference on
%        Hybrid Systems: Computation and Control.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork, nnHelper/validateNNoptions

% Authors:       Tobias Ladner
% Written:       28-March-2022
% Last update:   29-November-2022 (validateNNoptions)
%                16-February-2023 (re-organized structure)
%                21-February-2024 (moved internals to evaluate_)
% Last revision: 17-July-2023 (improved readability)

% ------------------------------ BEGIN CODE -------------------------------

% parse input
narginchk(2,5);

% validate parameters
[options, idxLayer] = setDefaultValues( ...
    {struct, 1:length(obj.layers)}, varargin);

% validate input
inputArgsCheck({ ...
    {obj, 'att', 'neuralNetwork'}; ...
    {input, 'att', {'numeric', 'interval', 'zonotope', 'polyZonotope', ...
        'taylm', 'conZonotope','gpuArray'}}; ...
    {options, 'att', 'struct'}; ...
    {idxLayer, 'att', 'numeric', 'vector'}; ...
    })
options = nnHelper.validateNNoptions(options);

% evaluate ----------------------------------------------------------------

r = evaluate_(obj,input,options,idxLayer);

end

% ------------------------------ END OF CODE ------------------------------
