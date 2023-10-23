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

evParams = struct;
evParams.refine = true;
evParams.reuse_bounds = true;

res = true;
for type = ["layer", "neuron"]
    for method = ["approx_error", "sensitivity", "both"]
        evParams.refine_type = type;
        evParams.refine_method = "approx_error";

        res = res && run_basic_nn_unittest(evParams);
    end
end

evParams.numRefine = 5;
res = res && run_basic_nn_unittest(evParams);

end

% ------------------------------ END OF CODE ------------------------------
