function res = test_nn_polyZonotope_quad()
% test_nn_polyZonotope_quad - tests nn using polyZonotopes with 'quad' 
%    approx
%
% Syntax:
%    res = test_nn_polyZonotope_quad()
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
options.nn.order = 2;

res = true;
for attention = ["ReLU", "sigmoid", "tanh"]
    options.nn.attention = attention;
    assert(res && run_basic_nn_unittest(options));
end

end

% ------------------------------ END OF CODE ------------------------------
