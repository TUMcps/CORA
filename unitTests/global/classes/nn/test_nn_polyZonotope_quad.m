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

evParams = struct;
evParams.order = 2;

res = true;
for attention = ["ReLU", "sigmoid", "tanh"]
    evParams.attention = attention;
    res = res && run_basic_nn_unittest(evParams);
end

end

% ------------------------------ END OF CODE ------------------------------
