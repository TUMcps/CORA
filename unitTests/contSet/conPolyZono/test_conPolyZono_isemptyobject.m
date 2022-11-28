function res = test_conPolyZono_isemptyobject
% test_conPolyZono_isemptyobject - unit test function of isemptyobject
%
% Syntax:  
%    res = test_conPolyZono_isemptyobject
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

% Author:       Mark Wetzlinger
% Written:      03-June-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% instantiate constrained polynomial zonotopes
cPZ1 = conPolyZono();

c = [0;0];
G = [1 0;0 1];
expMat = [1 0;0 1];
A = [1 -1];
b = 0;
expMat_ = [2 0;0 1];
cPZ2 = conPolyZono(c,G,expMat,A,b,expMat_);

% check results
res = isemptyobject(cPZ1) && ~isemptyobject(cPZ2);

%------------- END OF CODE --------------
