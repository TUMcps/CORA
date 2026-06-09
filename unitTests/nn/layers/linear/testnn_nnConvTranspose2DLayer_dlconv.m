function res = testnn_nnConvTranspose2DLayer_dlconv()
% testnn_nnConvTranspose2DLayer_dlconv - tests constructor and dlconv usage of nnConvTranspose2DLayer
%
% Syntax:
%    res = testnn_nnConvTranspose2DLayer_dlconv()
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

% Authors:       Benedikt Kellner
% Written:       01-December-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Setup options.
options.nn.use_dlconv = true;
options = nnHelper.validateNNoptions(options);

% Run test with deep learning tool box.
res = test_nn_nnConvTranspose2DLayer(options);

% ------------------------------ END OF CODE ------------------------------
