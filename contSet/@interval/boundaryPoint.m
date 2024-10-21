function x = boundaryPoint(I,dir,varargin)
% boundaryPoint - computes the point on the boundary of an interval along a
%    given direction, starting from a given start point, or, by default,
%    from the center of the set; for unbounded intervals, a start point
%    must be provided; any given start point must be contained in the set;
%    note that the vector may immediately reach the boundary of degenerate
%    intervals
%
% Syntax:
%    x = boundaryPoint(I,dir)
%    x = boundaryPoint(I,dir,startPoint)
%
% Inputs:
%    I - interval object
%    dir - direction
%    startPoint - start point for the direction vector
%
% Outputs:
%    x - point on the boundary of the interval
%
% Example: 
%    I = interval([1;-1],[4;3]);
%    dir = [1;1];
%    x = boundaryPoint(I,dir);
%    
%    figure; hold on;
%    plot(I);
%    plot(x(1),x(2),'.k');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/boundaryPoint

% Authors:       Mark Wetzlinger
% Written:       13-August-2024
% Last update:   17-October-2024 (MW, add start point)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check number of input arguments
narginchk(2,3);
% default starting point is the origin
startPoint = setDefaultValues({center(I)},varargin);
if any(isnan(startPoint))
    throw(CORAerror("CORA:wrongValue",'third',...
        'For unbounded sets, a start point must be provided'));
end

if all(withinTol(dir,0))
    throw(CORAerror("CORA:wrongValue",'second',...
        'Vector has to be non-zero.'));
end

% check for dimensions
equalDimCheck(I,dir);
equalDimCheck(I,startPoint);

% read out dimension
n = dim(I);

% for empty sets, return empty
if representsa_(I,'emptySet',0)
    x = zeros(n,0);
    return
end

% start point must be contained in the set
if ~contains(I,startPoint)
    throw(CORAerror('CORA:wrongValue','third',...
        'Start point must be contained in the set.'));
end

% shift interval by start point
I = I - startPoint;

% sign of direction
sign_dir = sign(dir);

% compute diameter
d = 2*rad(I);
d_inf = isinf(d);

% for positive/negative values, extract upper/lower bound
bound = infimum(I) + d .* max(sign_dir,zeros([n,1]));

% compute factor of limiting dimension
ratio = bound(dir ~= 0) ./ dir(dir ~= 0);

% take minimum value (excluding -Inf)
ratio = min(ratio(~isinf(ratio)));

% multiply (normalized) direction with that factor
x = dir * ratio + startPoint;

% dimensions with infinite width (unless direction is zero in that
% direction!)
inf_dirs = d_inf & ~withinTol(dir,0);
x(inf_dirs) = Inf * sign_dir(d_inf);

% ------------------------------ END OF CODE ------------------------------
