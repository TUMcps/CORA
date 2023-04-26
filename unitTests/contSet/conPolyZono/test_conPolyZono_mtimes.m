function res = test_conPolyZono_mtimes
% test_conPolyZono_mtimes - unit test function of mtimes
%
% Syntax:  
%    res = test_conPolyZono_mtimes
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

% Author:       Tobias Ladner
% Written:      06-April-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

resvec = true(0);

% init
c = [0;0];
G = [2 0 2; 0 2 2];
Grest = [1;6];
expMat = [1 0 1; 0 1 1; 0 0 0];
A = [2 2 4 -4];
b = 0;
expMat_ = [1 0 1 0; 0 1 1 0; 0 0 0 1];

% multiply with matrix
cPZ = conPolyZono(c,G,expMat,A,b,expMat_,Grest);
M = [1 2; 3 4];
McPZ = M * cPZ;
resvec(end+1) = compareMatrices(McPZ.c, M*c);
resvec(end+1) = compareMatrices(McPZ.G, M*G);
resvec(end+1) = compareMatrices(McPZ.Grest, M*Grest);
resvec(end+1) = compareMatrices(McPZ.expMat, expMat);
resvec(end+1) = compareMatrices(McPZ.A, A);
resvec(end+1) = compareMatrices(McPZ.b, b);
resvec(end+1) = compareMatrices(McPZ.expMat_, expMat_);

% test without Grest
cPZ = conPolyZono(c,G,expMat,A,b,expMat_);
M = [1 2; 3 4];
McPZ = M * cPZ;

resvec(end+1) = compareMatrices(McPZ.c, M*c);
resvec(end+1) = compareMatrices(McPZ.G, M*G);
resvec(end+1) = isempty(McPZ.Grest);

% test with only c
cPZ = conPolyZono(c);
M = [1 2; 3 4];
McPZ = M * cPZ;

resvec(end+1) = compareMatrices(McPZ.c, M*c);
resvec(end+1) = isempty(McPZ.G);
resvec(end+1) = isempty(McPZ.Grest);

% test empty case
resvec(end+1) = isemptyobject([] * conPolyZono());

% add results
res = all(resvec);

%------------- END OF CODE --------------