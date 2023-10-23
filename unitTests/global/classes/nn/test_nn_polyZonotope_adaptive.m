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

evParams = struct;

res = true;
for attention = ["ReLU", "sigmoid", "tanh"]
    evParams.attention = attention;
    res = res && run_basic_nn_unittest(evParams);
end

end

% ------------------------------ END OF CODE ------------------------------
