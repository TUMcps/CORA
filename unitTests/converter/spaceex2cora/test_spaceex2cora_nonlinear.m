function res = test_spaceex2cora_nonlinear
% test_spaceex2cora_nonlinear - test for model conversion from SpaceEx to
%    CORA for a simple nonlinear system without inputs
%
% Syntax:
%    test_spaceex2cora_nonlinear
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Mark Wetzlinger
% Written:       10-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;


% directory to SpaceEx model file
dir_spaceex = [CORAROOT filesep 'unitTests' filesep 'converter' ...
    filesep 'spaceex2cora' filesep 'testSystems'];

% convert SpaceEx model from .xml file
spaceex2cora([dir_spaceex filesep 'test_nonlinear1.xml']);
spaceex2cora([dir_spaceex filesep 'test_nonlinear2.xml']);

% instantiate system from converted SpaceEx model
sys_spaceex1 = test_nonlinear1;
sys_spaceex2 = test_nonlinear2;


% instantiate equivalent CORA models
f = @(x,u) [2*x(1)^2 + sqrt(x(2)); 4*exp(x(1)) - 3*sin(x(2))];
sys_cora1 = nonlinearSys('test_nonlinear1_St1_FlowEq',f);
f = @(x,u) [2*x(1)^2 + sqrt(x(2)) - cos(u(1)); 4*exp(x(1)) - 3*sin(x(2)) + log(u(2))];
sys_cora2 = nonlinearSys('test_nonlinear2_St1_FlowEq',f);


% compare systems
if sys_cora1 ~= sys_spaceex1
    res = false;
elseif sys_cora2 ~= sys_spaceex2
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
