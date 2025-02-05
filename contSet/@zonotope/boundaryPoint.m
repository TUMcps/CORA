function x = boundaryPoint(Z,dir,varargin)
% boundaryPoint - computes the point on the boundary of a zonotope along a
%    given direction, starting from a given start point, or, by default,
%    from the center of the set; note that the vector may immediately reach
%    the boundary of degenerate zonotopes
%
% Syntax:
%    x = boundaryPoint(Z,dir)
%    x = boundaryPoint(Z,dir,startPoint)
%
% Inputs:
%    Z - zonotope object
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
% See also: interval/boundaryPoint

% Authors:       Mark Wetzlinger
% Written:       13-April-2024
% Last update:   17-October-2024 (MW, algorithm for different start point)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check number of input arguments
narginchk(2,3);
% default starting point is the origin
startPoint = setDefaultValues({Z.c},varargin);
% ...center always exists, as zonotopes are bounded

if all(withinTol(dir,0))
    throw(CORAerror("CORA:wrongValue",'second',...
        'Vector has to be non-zero.'));
end

% check for dimensions
equalDimCheck(Z,dir);
equalDimCheck(Z,startPoint);

% read out dimension
n = dim(Z);

% empty set
if representsa_(Z,'emptySet',0)
    x = zeros(n,0);
    return
end

% start point must be contained in the set
if ~contains_(Z,startPoint,'exact',1e-8,0,false,false)
    throw(CORAerror('CORA:wrongValue','third',...
        'Start point must be contained in the set.'));
end

% 1D zonotope (start point has no effect)
if n == 1
    if dir < 0
        x = Z.c - sum(abs(Z.G));
    elseif dir > 0
        x = Z.c + sum(abs(Z.G));
    end
    return
end

% use zonotope norm if start point is the center (default case)
if all(withinTol(Z.c,startPoint,1e-12))
    % translate by the center
    Z_atorigin = zonotope(zeros(n,1),Z.G);
    
    % compute boundary point using the zonotope norm
    x = Z.c + dir ./ zonotopeNorm(Z_atorigin,dir);
    return
end


% if start point is another point contained in the zonotope, we formulate
% the boundary point computation as a linear program:
%   max t  s.t.  startPoint + t*dir \in Z
% which can be formulated as
%   min_{t,beta} -t 
%   s.t.  s + l*t = c + G beta  <=>  l*t - G beta = c - s
%         t >= 0, -1 <= beta <= 1
numGen = size(Z.G,2);

problem.f = [-1; zeros(numGen,1)];

problem.Aeq = [dir, -Z.G];
problem.beq = Z.c - startPoint;

% alternative to using lb/ub
% problem.Aineq = [-1, zeros(1,numGen); ...
%                  zeros(numGen,1), eye(numGen); ...
%                  zeros(numGen,1), -eye(numGen)];
% problem.bineq = [0; ones(numGen,1); -ones(numGen,1)];
problem.Aineq = [];
problem.bineq = [];

problem.lb = [0; -ones(numGen,1)];
problem.ub = [Inf; ones(numGen,1)];

% we should not run into any problems as the start point is guaranteed to
% be contained within the zonotope (by the check above)
[~,fval] = CORAlinprog(problem);

% compute solution (note the -1*fval because of optimizing -t)
x = startPoint - fval*dir;

% ------------------------------ END OF CODE ------------------------------
