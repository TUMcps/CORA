function res = isFullDim(P)
% isFullDim - checks if the dimension of the affine hull of a polytope is
%    equal to the dimension of its ambient space
%
% Syntax:  
%    res = isFullDim(P)
%
% Inputs:
%    P - mptPolytope object
%
% Outputs:
%    res - true/false
%
% Example:
%    P1 = mptPolytope([1 0;0 1;-1 0; 0 -1],[1;1;0;0]);
%    P2 = mptPolytope([1 0;0 1;-1 0; 0 -1],[1;0;0;0]);
%
%    isFullDim(P1)
%    isFullDim(P2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/isFullDim

% Author:       Niklas Kochdumper
% Written:      02-January-2020 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% check for emptiness
if isempty(P)
    res = false; return
end

% init result
res = true;

% get object properties
A = P.P.A;
b = P.P.b;

% normalize constraint matrix
n = sqrt(sum(A.^2,2));
A = A./n;
b = b./n;

% loop over all constraints
for i = 1:length(b)
    
    [res_,ind] = ismembertol(-A(i,:),A,eps,'ByRows',true);
    
    if res_ && b(ind) == -b(i)
        res = false; return;
    end
end

%------------- END OF CODE --------------