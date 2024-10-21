function res = test_nn_polyZonotope_lin()
% test_nn_polyZonotope_lin - tests nn using polyZonotopes with 'lin' 
%    approx
%
% Syntax:
%    res = test_nn_polyZonotope_lin()
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

% Authors:       Tobias Ladner
% Written:       24-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

options = struct;
options.nn.poly_method = "singh";

res = true;
for attention = ["ReLU", "sigmoid", "tanh"]
    options.nn.attention = attention;
    assert(res && run_basic_nn_unittest(options));
end

end

% ------------------------------ END OF CODE ------------------------------
