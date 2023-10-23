function res = test_nn_nnActivationLayer_getDerBounds()
% test_nn_nnActivationLayer_getDerBounds - tests if the computed bounds
%    are within the global derivative bounds
%
% Syntax:
%    res = test_nn_nnActivationLayer_getDerBounds()
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
% Written:       17-February-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

I = interval(0, 1);

layer = nnReLULayer();
res = res && I.contains(layer.getDerBounds(-1, 1));

layer = nnSigmoidLayer();
res = res && I.contains(layer.getDerBounds(-1, 1));

layer = nnTanhLayer();
res = res && I.contains(layer.getDerBounds(-1, 1));


% ------------------------------ END OF CODE ------------------------------
