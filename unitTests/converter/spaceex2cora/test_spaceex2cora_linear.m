function res = test_spaceex2cora_linear
% test_spaceex2cora_linear - test for model conversion from SpaceEx to CORA
%    for a simple linear system
%
% Syntax:
%    test_spaceex2cora_linear
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
spaceex2cora([dir_spaceex filesep 'test_linear1.xml']);
spaceex2cora([dir_spaceex filesep 'test_linear2.xml']);

% instantiate system from converted SpaceEx model
sys_spaceex1 = test_linear1;
sys_spaceex2 = test_linear2;


% instantiate equivalent CORA models
A = [2 1; 4 -3];
B = [0; 0];
sys_cora1 = linearSys(A,B);

A = [2 1; 4 -3];
B = [0.5; -0.25];
sys_cora2 = linearSys(A,B);

% compare systems
if sys_cora1 ~= sys_spaceex1
    res = false;
elseif sys_cora2 ~= sys_spaceex2
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
