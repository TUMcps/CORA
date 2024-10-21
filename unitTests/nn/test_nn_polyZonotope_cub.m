function res = test_nn_polyZonotope_cub()
% test_nn_polyZonotope_cub - tests nn using polyZonotopes with 'cub' approx
%    
%
% Syntax:
%    res = test_nn_polyZonotope_cub()
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
options.nn.order = 3;

res = true;
for attention = ["ReLU", "sigmoid", "tanh"]
    options.nn.attention = attention;
    assert(res && run_basic_nn_unittest(options));
end

end

% ------------------------------ END OF CODE ------------------------------
