function x = boundaryPoint(Z,dir)
% boundaryPoint - computes the point on the boundary of a zonotope along a
%    given direction; all-zero directions are not supported; the boundary
%    point of a degenerate zonotope along all directions that are not in
%    the affine hull of the zonotope is set to be the center
%
% Syntax:
%    x = boundaryPoint(Z,dir)
%
% Inputs:
%    Z - zonotope object
%    dir - direction
%
% Outputs:
%    x - point on the boundary of the zonotope
%
% Example: 
%    Z = zonotope([1;-1],[-3 2 1; -1 0 3]);
%    dir = [1;1];
%    x = boundaryPoint(Z,dir);
%    
%    figure; hold on;
%    plot(Z);
%    plot(x(1),x(2),'.k');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Mark Wetzlinger
% Written:       13-April-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check for dimensions
equalDimCheck(Z,dir);

if all(withinTol(dir,0))
    throw(CORAerror("CORA:notSupported","Vector has to be non-zero."));
end

% read out dimension
n = dim(Z);

% empty set
if representsa_(Z,'emptySet',0)
    x = zeros(n,0);
    return
end

% 1D zonotope
if n == 1
    if dir < 0
        x = Z.c - sum(abs(Z.G));
    elseif dir > 0
        x = Z.c + sum(abs(Z.G));
    end
    return
end

% translate by the center
Z_atorigin = zonotope(zeros(n,1),Z.G);

% compute boundary point using the zonotope norm
x = Z.c + dir ./ zonotopeNorm(Z_atorigin,dir);

% ------------------------------ END OF CODE ------------------------------
