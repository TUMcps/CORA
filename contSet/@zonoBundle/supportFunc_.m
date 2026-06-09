function [val,x] = supportFunc_(zB,dir,type,varargin)
% supportFunc_ - Calculate the upper or lower bound of a zonotope bundle
%    along a certain direction
%
% Syntax:
%    val = supportFunc_(zB,dir)
%    [val,x] = supportFunc_(zB,dir,type)
%
% Inputs:
%    zB - zonoBundle object
%    dir - direction for which the bounds are calculated (vector of size
%          (n,1) )
%    type - upper bound, lower bound, or both ('upper','lower','range')
%
% Outputs:
%    val - bound of the zonotope bundle in the specified direction
%    x - support vector
%
% Example: 
%    Z1 = zonotope([0 1 2 0;0 1 0 2]);
%    Z2 = zonotope([3 -0.5 3 0;-1 0.5 0 3]);
%    zB = zonoBundle({Z1,Z2});
%    val = supportFunc(zB,[1;1]);
%   
%    figure; hold on;
%    plot(Z1); plot(Z2); plot(zB);
%    plot(polytope([],[],[1,1],val),[1,2],'g');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/supportFunc, conZonotope/supportFunc_

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       19-November-2019
% Last update:   23-April-2023 (MW, fix empty case)
% Last revision: 27-March-2023 (MW, rename supportFunc_)

% ------------------------------ BEGIN CODE -------------------------------

% dimension
n = dim(zB);

% initialization
Aeq = []; beq = [];
lb = []; ub = [];

% loop over all parallel sets
for i = 1:zB.parallelSets
    
   % get object properties 
   Z = zB.Z{i};
   c = center(Z);
   G = generators(Z);
   nrGens = size(G,2);
   
   % construct equality constraint matrices
   Aeq = blkdiag(Aeq,-G);
   beq = [beq;c];
   lb = [lb;-ones(nrGens,1)];
   ub = [ub;ones(nrGens,1)];

end

% add optimal point as an additional variable
A = [eye(size(Aeq,2));-eye(size(Aeq,2))];
problem.Aineq = [zeros(size(A,1),n),A];
problem.bineq = [ub;-lb];

problem.Aeq = [repmat(eye(n),[zB.parallelSets,1]),Aeq];
problem.beq = beq;

f = [dir;zeros(length(lb),1)];

% solve LP via shared helper (handles lower/upper/range with x0 warm-start)
[val,x] = supportFunc_linprog(f, problem.Aineq, problem.bineq, ...
    problem.Aeq, problem.beq, [], [], type);

if nargout > 1
    % truncate support vector (remove generator factor variables)
    if strcmp(type,'range')
        if ~isempty(x)
            x = [x(1:n,1), x(1:n,2)];
        end
    else
        if ~isempty(x)
            x = x(1:n);
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
