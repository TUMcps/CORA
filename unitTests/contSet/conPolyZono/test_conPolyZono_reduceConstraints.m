function res = test_conPolyZono_reduceConstraints
% test_conPolyZono_reduceConstraints - unit test function for the
%    constraint reduction of constrained polynomial zonotopes
%
% Syntax:  
%    res = test_conPolyZono_reduceConstraints()
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
% See also: conPolyZono/reduceConstraints

% Author:       Niklas Kochdumper
% Written:      26-January-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

% Analytical Tests ----------------------------------------------------

% define constrained polynomial zonotope
c = [0;0];
G = [1 0 1 -1; 0 1 1 1];
expMat = [1 0 1 2; 0 1 1 0; 0 0 1 1];
A = [1 -0.5 0.5];
b = 0.5;
expMat_ = [0 1 2; 1 0 0; 0 1 0];

cPZ = conPolyZono(c,G,expMat,A,b,expMat_);

% remove constraints
cPZ = reduceConstraints(cPZ,0);

% compare with exact solution
c = [0;0.5];
G = [1  -1  0.5 0   0.5 -0.5;
     0   1  1  -0.5 0.5 -0.5];
expMat = [1 2 1 2 2 3;
          0 1 1 0 2 1];
id = [1;3];

cPZ_ = conPolyZono(c,G,expMat,[],[],[],[],id);

if ~isequal(cPZ,cPZ_)
    res = false;
end

%------------- END OF CODE --------------