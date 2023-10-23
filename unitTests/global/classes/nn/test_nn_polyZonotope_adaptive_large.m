function res = test_nn_polyZonotope_adaptive_large()
% test_nn_polyZonotope_adaptive_large - test large network
%
%
% Syntax:
%    res = test_nn_polyZonotope_adaptive_large
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
evParams.neurons = [2, 256, 256, 128, 64, 2];
evParams.sigmoid = "sigmoid";

res = run_basic_nn_unittest(evParams);

end

% ------------------------------ END OF CODE ------------------------------
