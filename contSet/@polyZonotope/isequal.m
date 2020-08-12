function res = isequal(pZ1,pZ2,tol)
% isequal - checks if two polyZonotopes are equal
%
% Syntax:  
%    res = isequal(Z1,Z2,tol)
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
%    pZ1 = polyZonotope(zeros(2,1),rand(3,5),[],randi([-2 2],3,5));
%    pZ2 = polyZonotope(zeros(2,1),rand(3,5),[],randi([-2 2],3,5));
%    isequal(pZ1,pZ2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Mark Wetzlinger
% Written:       01-May-2020
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

if nargin == 2
    tol = eps;
end

res = false;

% compare dimensions (quick check)
if dim(pZ1) ~= dim(pZ2)
    return
end

% compare centers (quick check)
if any(abs(center(pZ1) - center(pZ2)) > tol)
    return
end

% remove redundant exponents and delete zero-length generators
expMat1 = pZ1.expMat;
G1 = pZ1.G;
[expMat1new,G1new] = removeRedundantExponents(expMat1,G1);
pZ1 = deleteZeros(polyZonotope(pZ1.c,G1new,pZ1.Grest,expMat1new));
expMat2 = pZ2.expMat;
G2 = pZ2.G;
[expMat2new,G2new] = removeRedundantExponents(expMat2,G2);
pZ2 = deleteZeros(polyZonotope(pZ2.c,G2new,pZ2.Grest,expMat2new));

% compare sizes of exponent matrices
if any(size(pZ1.expMat) - size(pZ2.expMat))
    return
end

% compare dependent generators
if ~(all(all(abs(pZ1.G - pZ2.G) < tol)))
    return
end

% compare independent generators
if ~(all(all(abs(pZ1.Grest - pZ2.Grest) < tol)))
    return
end

res = true;

%------------- END OF CODE --------------

