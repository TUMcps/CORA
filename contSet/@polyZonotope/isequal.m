function res = isequal(pZ1,pZ2,varargin)
% isequal - checks if two polynomial zonotopes are equal
%
% Syntax:  
%    res = isequal(pZ1,pZ2)
%    res = isequal(pZ1,pZ2,tol)
%
% Inputs:
%    pZ1 - polyZonotope object
%    pZ2 - polyZonotope object
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
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

% too many input arguments
if nargin > 3
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

% parse input arguments
tol = setDefaultValues({eps},varargin);

% check input arguments
inputArgsCheck({{pZ1,'att','polyZonotope'};
                {pZ2,'att','polyZonotope'};
                {tol,'att','numeric',{'nonnan','scalar','nonnegative'}}});

% assume false
res = false;

% remove redundancies in representation
pZ1 = compact(pZ1);
pZ2 = compact(pZ2);

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

% jointly compare dependent generators and exponent matrices
if ~compareMatrices([pZ1.G;E1],[pZ2.G;E2], tol)
    return
end
if ~compareMatrices(pZ1.Grest,pZ2.Grest, tol)
    return
end
if ~all(withinTol(pZ1.c, pZ2.c, tol))
    return
end

% all checks ok
res = true;

%------------- END OF CODE --------------