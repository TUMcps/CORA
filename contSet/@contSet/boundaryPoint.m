function x = boundaryPoint(S,dir,varargin)
% boundaryPoint - computes the point on the boundary of a set along a
%    given direction, starting from a given start point, or, by default,
%    from the center of the set; for unbounded sets, a start point must be
%    provided; any given start point must be contained in the set; note
%    that the vector may immediately reach the boundary of degenerate sets
%
% Syntax:
%    x = boundaryPoint(S,dir)
%    x = boundaryPoint(S,dir,startPoint)
%
% Inputs:
%    S - contSet object
%    dir - direction
%    startPoint - start point for the direction vector
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
% See also: none

% Authors:       Mark Wetzlinger
% Written:       17-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check number of input arguments
narginchk(2,3);
% default starting point is the center: caution, as this may be
% time-consuming to compute
startPoint = setDefaultValues({center(S)},varargin);
if any(isnan(startPoint))
    throw(CORAerror("CORA:wrongValue",'third',...
        'For unbounded sets, a start point must be provided'));
end

if all(withinTol(dir,0))
    throw(CORAerror("CORA:wrongValue",'second',...
        'Vector has to be non-zero.'));
end

% check for dimensions
equalDimCheck(S,dir);
equalDimCheck(S,startPoint);

% for empty sets, return empty
if representsa_(S,'emptySet',0)
    x = zeros(dim(S),0);
    return
end

% start point must be contained in the set
if ~contains(S,startPoint)
    throw(CORAerror('CORA:wrongValue','third',...
        'Start point must be contained in the set.'));
end

% not implemented in subclasses
throw(CORAerror('CORA:noops',S));

% ------------------------------ END OF CODE ------------------------------
