function res = test_cora2spaceex_nonlinear
% test_cora2spaceex_nonlinear - test for model conversion from CORA to 
%   SpaceEx for a nonlinear system
%
% Syntax:
%    test_cora2spaceex_nonlinear
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

% load model
name = 'conv_test_cora2spaceex_nonlinear_St1_FlowEq';
sys = nonlinearSys(name,@doublePendulum);

% convert model to SpaceEx format
cora2spaceex(sys,'model_test_cora2spaceex_nonlinear');

% import model from SpaceEx format
spaceex2cora('model_test_cora2spaceex_nonlinear',[],[],...
                                    'conv_test_cora2spaceex_nonlinear');
sys_ = conv_test_cora2spaceex_nonlinear();

% compare hybrid automata
assert(sys == sys_);

% ------------------------------ END OF CODE ------------------------------
