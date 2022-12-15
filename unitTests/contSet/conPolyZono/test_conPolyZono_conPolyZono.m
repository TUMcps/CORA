function res = test_conPolyZono_conPolyZono
% test_conPolyZono_conPolyZono - unit test function of constructor
%
% Syntax:  
%    res = test_conPolyZono_conPolyZono
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

% Author:       Tobias Ladner
% Written:      15-December-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% test example in doctring
c = [0;0];
G = [1 0 1 -1; 0 1 1 1];
expMat = [1 0 1 2; 0 1 1 0; 0 0 1 1];
A = [1 -0.5 0.5];
b = 0.5;
expMat_ = [0 1 2; 1 0 0; 0 1 0];

cPZ = conPolyZono(c,G,expMat,A,b,expMat_);

% test variants in syntax of docstring
Grest = [4 1; 0 2];
id = [1;2;4];

cPZ = conPolyZono(c,G,expMat);
cPZ = conPolyZono(c,G,expMat,Grest);
cPZ = conPolyZono(c,G,expMat,Grest,id);
cPZ = conPolyZono(c,G,expMat,A,b,expMat_);
cPZ = conPolyZono(c,G,expMat,A,b,expMat_,Grest);
cPZ = conPolyZono(c,G,expMat,A,b,expMat_,Grest,id);

% test copy constructor
cPZ = conPolyZono(cPZ);

% test empty set
cPZ = conPolyZono();

res = true;

%------------- END OF CODE --------------
