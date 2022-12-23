function res = isequal(cZ,S,varargin)
% isequal - checks if a constrained zonotope represents the same set as
%    another set
%
% Syntax:  
%    res = isequal(cZ,S)
%    res = isequal(cZ,S,tol)
%
% Inputs:
%    cZ - conZonotope object
%    S - contSet object
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example:
%    c = [0;0]; G = [1 0 1;0 1 1];
%    A = [1 1 1]; b = 4;
%    cZ = conZonotope([c,G],A,b);
%
%    isequal(cZ,cZ + [0;eps],0)
%    isequal(cZ,cZ + [0;eps],1e-14)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Mark Wetzlinger
% Written:      19-December-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% too many input arguments
if nargin > 3
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

% set default tolerance
tol = setDefaultValues({1e-12},varargin);

% check input arguments
inputArgsCheck({ {cZ,'att','conZonotope'}, ...
                 {S,'att','contSet'}, ...
                 {tol,'att','numeric',{'nonempty','scalar','finite','nonnegative'}}});

% only implemented for conZonotope - conZonotope/zonotope
if ~( isa(S,'conZonotope') || isa(S,'zonotope') )
    throw(CORAerror('CORA:noops',cZ,S));
end

% second set: zonotope
if isa(S,'zonotope')
    % first set: conZonotope without constraints
    if isempty(cZ.A) && isempty(cZ.b)
        res = isequal(zonotope(cZ),S,tol);
    else
        % constraints given
        res = false;
    end
    return
end

% two constrained zonotopes

% centers have to be the same
if ~withinTol(cZ.Z(:,1),S.Z(:,1),tol)
    res = false;
    return
end

% check if both are represented completely equally
if compareMatrices([cZ.Z(:,2:end); cZ.A],[S.Z(:,2:end); S.A],tol) ...
        && compareMatrices([cZ.A cZ.b],[S.A S.b],tol)
    res = true;
    return
end

% last resort: check vertices (computationally expensive!)
res = compareMatrices(vertices(cZ),vertices(S),tol);

%------------- END OF CODE --------------