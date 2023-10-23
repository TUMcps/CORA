function S = compact(S,varargin)
% compact - removes redundancies in the representation of a set, the
%    resulting set is equal to the original set, but minimal in its
%    representation
%
% Syntax:
%    S = compact(S)
%    S = compact(S,method)
%    S = compact(S,method,tol)
%
% Inputs:
%    S - contSet object
%    method - (optional) method for redundancy removal
%    tol - (optional) tolerance
%
% Outputs:
%    S - contSet object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       29-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% no changes to following classes, which are always in their minimal
% representation:
% - capsule, ellipsoid, emptySet, fullspace, halfspace, interval
if isa(S,'capsule') || isa(S,'ellipsoid') || isa(S,'emptySet') ...
        || isa(S,'fullspace') || isa(S,'halfspace') || isa(S,'interval')
    return;
end

% parse input arguments, set default values
if isa(S,'zonotope')
    [method,tol] = setDefaultValues({'zeros',eps},varargin);
    methods = {'all','zeros','aligned'};
    % reset tolerance for 'aligned' method
    if nargin < 3 && strcmp(method,'aligned')
        tol = 1e-3;
    end

elseif isa(S,'polytope')
    [method,tol] = setDefaultValues({'all',1e-9},varargin);
    methods = {'all','zeros','A','Ae','aligned','V'};

elseif isa(S,'conZonotope')
    [method,tol] = setDefaultValues({'all',eps},varargin);
    methods = {'all','zeros'};

elseif isa(S,'polyZonotope')
    [method,tol] = setDefaultValues({'all',eps},varargin);
    methods = {'all','states','exponentMatrix'};

elseif isa(S,'conPolyZono')
    [method,tol] = setDefaultValues({'all',eps},varargin);
    methods = {'all','states','constraints','exponentMatrix'};

elseif isa(S,'levelSet')
    [method,tol] = setDefaultValues({'all',eps},varargin);
    methods = 'all';

end
% check input arguments
inputArgsCheck({{method,'str',methods},...
    {tol,'att','numeric',{'scalar','nonnegative'}}});


% other classes may implement compact
try
    S = compact_(S,method,tol);
catch ME
    if strcmp(ME.identifier,'')
        throw(CORAerror('CORA:noops',S));
    else
        rethrow(ME);
    end
end

% ------------------------------ END OF CODE ------------------------------
