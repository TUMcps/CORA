function [val,x,ksi] = supportFunc_(cZ,dir,type,varargin)
% supportFunc_ - Calculate the upper or lower bound of a constrained 
%    zonotope along a certain direction
%
% Syntax:
%    val = supportFunc_(cZ,dir)
%    [val,x,ksi] = supportFunc_(cZ,dir)
%    [val,x,ksi] = supportFunc_(cZ,dir,type)
%
% Inputs:
%    cZ - conZonotope object
%    dir - direction for which the bounds are calculated (vector of size
%          (n,1) )
%    type - upper bound, lower bound, or both ('upper','lower','range')
%
% Outputs:
%    val - bound of the constrained zonotope in the specified direction
%    x - support vector (column vectors if type = 'range')
%    ksi - factor values that correspond to the bound (column vectors if
%          type = 'range')
%
% Example:
%    Z = [0 1 0 1;0 1 2 -1];
%    A = [-2 1 -1]; b = 2;
%    cZ = conZonotope(Z,A,b);
%
%    dir = [1;1];
%    val = supportFunc(cZ,dir);
%
%    figure; hold on; xlim([-3,1]); ylim([-3,4]);
%    plot(cZ);
%    plot(polytope([],[],dir',val),[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/supportFunc, zonotope/supportFunc_

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       22-May-2018
% Last update:   10-December-2022 (MW, add type = 'range')
%                16-December-2022 (MW, simplify code, add call to MOSEK)
% Last revision: 27-March-2023 (MW, rename supportFunc_)

% ------------------------------ BEGIN CODE -------------------------------

% non-constrained case
if isempty(cZ.A) || all(all(cZ.A == 0))
    [val,x,ksi] = aux_supportFunc_unconstrained(cZ,dir,type);
    return
end

% project the zonotope onto the direction
G_proj = dir' * cZ.G;

% check if one of the special cases are applicable
[fval,x,ksi] = aux_trySpecialCases(G_proj,cZ.A,cZ.b);

empty = false;
if isempty(fval)
    % linear program:
    % max_{x \in cZ} dir' * x
    % s.t. x = c + G*ksi
    %      A * ksi = b
    [ksi,fval,val,empty] = aux_evaluateLP(G_proj,cZ.A,cZ.b,type);
end

% unless emptiness determined, calculate bound by adding the zonotope center
if ~empty
    if strcmp(type,'range')
        temp = dir' * cZ.c + fval;
        val = interval(temp(1),temp(2));
    elseif any(strcmp(type,{'upper','lower'}))
        val = dir' * cZ.c + fval;
    end
    
    % calculate support vector
    if nargout >= 2
        x = cZ.c + cZ.G*ksi;
    end
end

end


% Auxiliary functions -----------------------------------------------------

function [ksi,fval,val,empty] = aux_evaluateLP(G_proj,A,b,type)

% number of generators -> number of optimization variables
nrGen = numel(G_proj);

problem.Aineq = [];
problem.bineq = [];
problem.Aeq = A;
problem.beq = b;
% ksi in [-1, 1]
problem.lb = -ones(nrGen,1);
problem.ub = ones(nrGen,1);

% for easier concatenation if type == 'range'
fval = []; ksi = []; val = []; empty = false;

% lower bound
if any(strcmp(type,{'lower','range'}))
    problem.f = G_proj';
    [ksi,fval,exitflag] = CORAlinprog(problem);

    if exitflag == -2
        % primal infeasible -> empty set
        empty = true;
        if strcmp(type,'range')
            val = interval.empty(1);
        elseif strcmp(type,'lower')
            val = Inf;
        end
        return
    elseif exitflag <= -3
        % should not be unbounded, or other solver issue...
        throw(CORAerror('CORA:solverIssue'));
    end
end

% upper bound: since linprog always computes a minimization, we
% need to multiply with -1 in some parts...
if any(strcmp(type,{'upper','range'}))
    problem.f = -G_proj';
    [ksi_,fval_,exitflag] = CORAlinprog(problem);

    if exitflag == -2
        % primal infeasible -> empty set
        empty = true;
        if strcmp(type,'range')
            val = interval.empty(1);
        elseif strcmp(type,'upper')
            val = -Inf;
        end
        return
    elseif exitflag == 1
        % combine factors
        ksi = [ksi, ksi_];
        fval = [fval, -fval_];
    elseif exitflag <= -3
        % should not be unbounded, or other solver issue...
        throw(CORAerror('CORA:solverIssue'));
    end
end

end

function [val,x,ksi] = aux_supportFunc_unconstrained(cZ,dir,type)

% project zonotope onto the given direction
cZ_proj = dir'*cZ;
I = interval(cZ_proj);

% trivially empty case
if representsa_(I,'emptySet',0)
    x = []; ksi = [];
    if strcmp(type,'upper')
        val = -Inf;
    elseif strcmp(type,'lower')
        val = Inf;
    elseif strcmp(type,'range')
        val = interval(-Inf,Inf);
    end
    return
end

% determine bound(s) and support vector
if strcmp(type,'upper')
    val = supremum(I);
    ksi = sign(cZ_proj.G)';
    x = cZ.c + cZ.G*ksi;
elseif strcmp(type,'lower')
    val = infimum(I); 
    ksi = -sign(cZ_proj.G)';
    x = cZ.c + cZ.G*ksi;
elseif strcmp(type,'range')
    val = I;
    ksi = [-sign(cZ_proj.G)', sign(cZ_proj.G)'];
    x = [];
end

end

function [fval,x,ksi] = aux_trySpecialCases(G_proj,A,b)

fval = [];
x = [];
ksi = [];

nrGen = numel(G_proj);

% special cases
if ~any(G_proj)
    % projection of the zonotope equals 0 -> no extension of the set in
    % this direction
    fval = 0;
    ksi = zeros(nrGen,1);

elseif nargout == 1

    % check if the objective function is identical to the normal vector of
    % a constraint -> only one solution
    idxPos = all(withinTol(G_proj,A),2);
    idxNeg = all(withinTol(-G_proj,A),2);
    % if more than one constraint eligible, skip this variant and go for
    % linear program evaluation
    if nnz(idxPos) + nnz(idxNeg) == 1
        if any(idxPos)
            fval = b(idxPos);
        elseif any(idxNeg)
            fval = -b(idxNeg);
        end
        if strcmp(type,'range')
            fval = [fval fval];
        end
    end
end

end

% ------------------------------ END OF CODE ------------------------------
