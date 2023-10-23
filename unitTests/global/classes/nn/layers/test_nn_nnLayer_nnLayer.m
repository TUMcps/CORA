function [res] = test_nn_nnLayer_nnLayer()
% test_nn_nnLayer_nnLayer - tests constructor of nnLayer
%
% Syntax:
%    res = test_nn_nnLayer_nnLayer()
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
% Written:       23-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% check construction of layers
layer = nnSigmoidLayer();
layer = nnTanhLayer();
layer = nnReLULayer();

% check if name is set correctly
name = "TestLayer";
layer = nnSigmoidLayer(name);
layer = nnTanhLayer(name);
layer = nnReLULayer(name);

if layer.name ~= name
    res = false;
end

end

% ------------------------------ END OF CODE ------------------------------
