function cZ = project(cZ,dims)
% project - projects a constrained zonotope onto the specified dimensions
%
% Syntax:
%    cZ = project(cZ,dims)
%
% Inputs:
%    cZ - (conZonotope) constrained zonotope
%    dims - dimensions for projection
%
% Outputs:
%    cZ - (conZonotope) projected constrained zonotope
%
% Example: 
%    Z = [0 1 0 1 0;0 1 2 -1 0;0 0 0 0 1];
%    A = [-2 1 -1 0]; b = 2;
%    cZ = conZonotope(Z,A,b);
%    cZres = project(cZ,[1,2]);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       11-May-2018
% Last update:   14-March-2021 (MW, input args handling, empty set)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

try
    % try projection first, as isempty() can be quite costly
    cZ.c = cZ.c(dims,:);
    cZ.G = cZ.G(dims,:);
catch ME
    % now check if conZonotope was empty
    if representsa_(cZ,'emptySet',eps)
        throw(CORAerror('CORA:emptySet'));
    elseif any(dims > dim(cZ))
        throw(CORAerror('CORA:wrongValue','second','should not exceed dimension of conZonotope'));
    else
        rethrow(ME);
    end
end

% ------------------------------ END OF CODE ------------------------------
