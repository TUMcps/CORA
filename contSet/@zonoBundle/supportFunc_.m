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

problem.lb = [];
problem.ub = [];

% upper/lower bound or range
switch type 
    case 'lower'
    
        % solve linear program
        problem.f = f';
        [x,val,exitflag] = CORAlinprog(problem);
        if exitflag == -2
            % primal infeasible -> empty set
            val = Inf; x = []; return
        elseif exitflag <= -3
            % should not be unbounded, or other solver issue...
            throw(CORAerror('CORA:solverIssue'));
        end
    
    case 'upper'
    
        % solve linear program
        problem.f = -f';
        [x,val,exitflag] = CORAlinprog(problem);
        if exitflag == -2
            % primal infeasible -> empty set
            val = -Inf; x = []; return
        elseif exitflag <= -3
            % should not be unbounded, or other solver issue...
            throw(CORAerror('CORA:solverIssue'));
        end
        val = -val;
   
    case 'range'

        % solve linear program for upper bound
        problem.f = -f';
        [x_upper,val_upper,exitflag] = CORAlinprog(problem);
        if exitflag == -2
            % primal infeasible -> empty set
            val = interval(-Inf,Inf); x = []; return
        elseif exitflag <= -3
            % should not be unbounded, or other solver issue...
            throw(CORAerror('CORA:solverIssue'));
        end
        val_upper = -val_upper;
    
        % solve linear program for lower bound
        problem.f = f';
        [x_lower,val_lower] = CORAlinprog(problem);
    
        % combine results for output args
        val = interval(val_lower,val_upper);

end


if nargout > 1
    % truncate support vector
    if strcmp(type,'range')
        x = [x_lower(1:n), x_upper(1:n)];
    else
        x = x(1:n);
    end
end

% ------------------------------ END OF CODE ------------------------------
