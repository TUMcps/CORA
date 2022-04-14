function res = isequal(pZ1,pZ2,tol)
% isequal - checks if two polyZonotopes are equal
%
% Syntax:  
%    res = isequal(pZ1,pZ2)
%    res = isequal(pZ1,pZ2,tol)
%
% Inputs:
%    pZ1 - polyZonotope object
%    pZ2 - polyZonotope object
%    tol - tolerance (optional)
%
% Outputs:
%    res - boolean whether pZ1 and pZ2 are equal
%
% Example: 
%    pZ1 = polyZonotope([0;0],[1 0 1;0 -1 1],[0.4 0;0.1 1],[1 0 2;0 1 1]);
%    pZ2 = polyZonotope([0;0],[1 1 0;1 0 -1],[0 0.4;1 0.1],[2 1 0;1 0 1]);
%    isequal(pZ1,pZ2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/isequal

% Author:        Mark Wetzlinger
% Written:       01-May-2020
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

if nargin == 2
    tol = eps;
end

res = false;
pZ1 = removeRedundancies(pZ1);
pZ2 = removeRedundancies(pZ2);
% compare dimensions (quick check)
if dim(pZ1) ~= dim(pZ2)
    return
end

% remove redundancies
pZ1 = compact(pZ1);
pZ2 = compact(pZ2);

% compare number of generators (quick check)
if size(pZ1.G,2) ~= size(pZ2.G,2) || size(pZ1.Grest,2) ~= size(pZ2.Grest,2)
   return 
end

% compare identifier vectors
temp1 = sort(pZ1.id); temp2 = sort(unique([pZ1.id;pZ2.id]));
if length(temp1) ~= length(temp2) || ~all(temp1 == temp2)
   return;
elseif ~all(pZ1.id == pZ2.id)
   [~,E1,E2] = mergeExpMatrix(pZ1.id,pZ2.id,pZ1.expMat,pZ2.expMat);
else
   E1 = pZ1.expMat; E2 = pZ2.expMat;
end

% compare exponent matrices
if ~(all(all(abs(E1 - E2) < tol)))
    return
end

% compare dependent generators
if ~(all(all(abs(pZ1.G - pZ2.G) < tol)))
    return
end

% compare center and independent generators
Z1 = zonotope(pZ1.c,pZ1.Grest);
Z2 = zonotope(pZ2.c,pZ2.Grest);
res = isequal(Z1,Z2);

%------------- END OF CODE --------------