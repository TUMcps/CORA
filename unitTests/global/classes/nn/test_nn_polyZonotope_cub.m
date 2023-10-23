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

evParams = struct;
evParams.order = 3;

res = true;
for attention = ["ReLU", "sigmoid", "tanh"]
    evParams.attention = attention;
    res = res && run_basic_nn_unittest(evParams);
end

end

% ------------------------------ END OF CODE ------------------------------
