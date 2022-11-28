function res = isequal(Z1,Z2,varargin)
% isequal - checks if two zonotopes are equal (note: no deletion of aligned
%    generators since this is quite costly)
%
% Syntax:  
%    res = isequal(Z1,Z2,tol)
%
% Inputs:
%    Z1 - zonotope object
%    Z2 - zonotope object
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    Z1 = zonotope(zeros(2,1),[1 0 2; 2 0 1]);
%    Z2 = zonotope(zeros(2,1),[1 2 0; 2 1 0]);
%    isequal(Z1,Z2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      16-September-2019
% Last update:  09-June-2020 (include re-ordering of generators)
%               13-November-2022 (MW, integrate modular check)
% Last revision:---

%------------- BEGIN CODE --------------

% parse input arguments
tol = setDefaultValues({eps},varargin{:});

% check input arguments
inputArgsCheck({{Z1,'att','zonotope'};
                {Z2,'att','zonotope'};
                {tol,'att','numeric',{'nonnan','scalar','nonnegative'}}});

% init result
res = false;

% compare dimensions (quick check)
if dim(Z1) ~= dim(Z2)
    return
end

% compare centers (quick check)
if ~all(withinTol(center(Z1),center(Z2),tol))
    return
end

% delete zeros from generator matrices
G1 = generators(deleteZeros(Z1));
G2 = generators(deleteZeros(Z2));

% compare number of generators
if size(G1,2) ~= size(G2,2)
    return
end

% compare generator matrices
res = compareMatrices(G1,G2,tol);

%------------- END OF CODE --------------