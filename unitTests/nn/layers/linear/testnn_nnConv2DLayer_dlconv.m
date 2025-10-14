function res = testnn_nnConv2DLayer_dlconv()
% testnn_nnConv2DLayer_dlconv - tests constructor of nnConv2DLayer
%
% Syntax:
%    res = testnn_nnConv2DLayer_dlconv()
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
% See also: none

% Authors:       Lukas Koller
% Written:       10-April-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Setup options.
options.nn.use_dlconv = true;
options = nnHelper.validateNNoptions(options);

% Run test with deep learning tool box.
res = test_nn_nnConv2DLayer(options);

% ------------------------------ END OF CODE ------------------------------
