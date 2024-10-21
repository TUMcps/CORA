function res = test_nn_polyZonotope_adaptive_refine()
% test_nn_polyZonotope_adaptive_refine - test refinement
%
%
% Syntax:
%    res = test_nn_polyZonotope_adaptive_refine
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       24-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

options = struct;
options.nn.refine = true;
options.nn.reuse_bounds = true;

res = true;
for type = ["layer", "neuron"]
    for method = ["approx_error", "sensitivity", "both"]
        options.nn.refine_type = type;
        options.nn.refine_method = "approx_error";

        assert(res && run_basic_nn_unittest(options));
    end
end

options.numRefine = 5;
assert(res && run_basic_nn_unittest(options));

end

% ------------------------------ END OF CODE ------------------------------
