function res = isequal(P1,P2,varargin)
% isequal - check ifs two polytopes are equal
%
% Syntax:
%    res = isequal(P1,P2)
%    res = isequal(P1,P2,tol)
%
% Inputs:
%    P1 - polytope object 
%    P2 - polytope object 
%    tol - (optional) tolerance
%
% Outputs:
%    res - result of comparison
%
% Example: 
%    P1 = polytope([-1 -1; 1 0;-1 0; 0 1; 0 -1],[2;3;2;3;2]);
%    P2 = polytope([0 -1; 0 1;-1 0; 1 0; -1 -1],[2;3;2;3;2]);
%    res = isequal(P1,P2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       27-April-2023 (MW, migrated code from eq)
% Last update:   27-July-2023 (MW, handle 1D case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set default values
tol = setDefaultValues(1e-9,varargin);

%check dimensions
equalDimCheck(P1,P2);

% one-dimensional cases: compare vertices (fast) and compare
if dim(P1) == 1
    V1 = vertices_(P1,'lcon2vert');
    V2 = vertices_(P2,'lcon2vert');
    if xor(any(isinf(V1)),any(isinf(V2)))
        res = false;
    elseif size(V1,2) ~= size(V2,2) 
        res = false;
    elseif any(isinf(V1))
        res = all(V1 == V2);
    else
        res = compareMatrices(V1,V2,tol);
    end
    return
end

% check if P1 and P2 contain each other
res = contains_(P1,P2,'exact',tol) && contains_(P2,P1,'exact',tol);

% ------------------------------ END OF CODE ------------------------------
