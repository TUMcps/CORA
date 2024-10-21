function res = test_nn_polyZonotope_adaptive()
% test_nn_polyZonotope_adaptive - test 'adaptive'
%
%
% Syntax:
%    res = test_nn_polyZonotope_adaptive
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

res = true;
for attention = ["ReLU", "sigmoid", "tanh"]
    options.nn.attention = attention;
    assert(res && run_basic_nn_unittest(options));
end

end

% ------------------------------ END OF CODE ------------------------------
