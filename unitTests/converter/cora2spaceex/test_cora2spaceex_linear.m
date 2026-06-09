function res = test_cora2spaceex_linear
% test_cora2spaceex_linear - test for model conversion from CORA to 
%   SpaceEx for a linear system
%
% Syntax:
%    test_cora2spaceex_linear
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Niklas Kochdumper
% Written:       14-November-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
% assume true
res = true;

% generate random model
sys = linearSys(rand(2),rand(2,3),rand(2,1));

% convert model to SpaceEx format
cora2spaceex(sys,'model_test_cora2spaceex_linear');

% import model from SpaceEx format
spaceex2cora('model_test_cora2spaceex_linear',[],[],...
                                    'conv_test_cora2spaceex_linear');
sys_ = conv_test_cora2spaceex_linear();

% compare hybrid automata
assert(sys == sys_);

% ------------------------------ END OF CODE ------------------------------
