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
%    plot(conHyperplane(halfspace([1;1],val)),[1,2],'g');
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

% fully-empty zonoBundle
if representsa_(zB,'emptySet',1e-8)
    x = [];
    if strcmp(type,'upper')
        val = -Inf;
    elseif strcmp(type,'lower')
        val = Inf;
    elseif strcmp(type,'range')
        val = interval(-Inf,Inf);
    end
    return
end

% initialization
Aeq = [];
beq = [];
lb = [];
ub = [];

% loop over all parallel sets
for i = 1:zB.parallelSets
    
   % get object properties 
   Z = zB.Z{i};
   G = generators(Z);
   c = center(Z);
   [n,nrGens] = size(G);
   
   % construct equality constraint matrices
   Aeq = blkdiag(Aeq,-G);
   beq = [beq;c];
   lb = [lb;-ones(nrGens,1)];
   ub = [ub;ones(nrGens,1)];

end

% add optimal point as an additional variable
A = [eye(size(Aeq,2));-eye(size(Aeq,2))];
A = [zeros(size(A,1),n),A];

b = [ub;-lb];

Aeq = [repmat(eye(n),[zB.parallelSets,1]),Aeq];

f = [dir;zeros(length(lb),1)];

% linear program options
persistent options
if isempty(options)
    options = optimoptions('linprog','display','off');
end

% upper or lower bound
if strcmp(type,'lower')
    
    % solve linear program
    [x,val,exitflag] = linprog(f',A,b,Aeq,beq,[],[],options);
    if exitflag == -2
        % primal infeasible -> empty set
        val = Inf; x = []; return
    elseif exitflag <= -3
        % should not be unbounded, or other solver issue...
        throw(CORAerror('CORA:solverIssue'));
    end
    
elseif strcmp(type,'upper')
    
    % solve linear program
    [x,val,exitflag] = linprog(-f',A,b,Aeq,beq,[],[],options);
    if exitflag == -2
        % primal infeasible -> empty set
        val = -Inf; x = []; return
    elseif exitflag <= -3
        % should not be unbounded, or other solver issue...
        throw(CORAerror('CORA:solverIssue'));
    end
    val = -val;
   
elseif strcmp(type,'range')

    % solve linear program for upper bound
    [x_upper,val_upper,exitflag] = linprog(-f',A,b,Aeq,beq,[],[],options);
    if exitflag == -2
        % primal infeasible -> empty set
        val = interval(-Inf,Inf); x = []; return
    elseif exitflag <= -3
        % should not be unbounded, or other solver issue...
        throw(CORAerror('CORA:solverIssue'));
    end
    val_upper = -val_upper;
    % solve linear program for lower bound
    [x_lower,val_lower] = linprog(f',A,b,Aeq,beq,[],[],options);

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
