function res = test_polyZonotope_mtimes
% test_polyZonotope_mtimes - unit test function for multiplication between
%    an interval matrix and a zonotope 
%
% Syntax:  
%    res = test_polyZonotope_mtimes
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

% Author:       Niklas Kochdumper
% Written:      26-June-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% create polynomial zonotope
c = [1;2];
G = [1 2 1 -3; 1 -1 2 -1];
expMat = [1 0 0 2; 0 1 2 1];
Grest = [];
pZ = polyZonotope(c,G,Grest,expMat);

% create interval matrix
matrix = interval([1 -0.5; -1 0], [3 0.5; 3 2]);

% multiply interval matrix with polynomial zonotope
pZres = matrix * pZ;

% define ground truth
c = [2; 3];
G =  [2 4 2 -6; 2 1 3 -4];
Grest = [11.5 0; 0 23];

% check for correctness
if ~all(withinTol(c,pZres.c)) || ~compareMatrices(G,pZres.G) ...
        || ~compareMatrices(Grest,pZres.Grest)
    throw(CORAerror('CORA:testFailed'));
end

% empty set
% if ~isempty(matrix*polyZonotope())
%     throw(CORAerror('CORA:testFailed'));
% end

res = true;

%------------- END OF CODE --------------