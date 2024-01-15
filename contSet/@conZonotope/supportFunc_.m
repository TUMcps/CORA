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
%    plot(conHyperplane(dir,val),[1,2],'r');
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

% empty case
if representsa_(cZ,'emptySet',0)
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

% non-constrained case
if isempty(cZ.A) || all(all(cZ.A == 0)) 
    
    % project zonotope onto the direction
    temp = dir'*cZ;
    I = interval(temp);
    
    % determine upper or lower bound
    if strcmp(type,'upper')
        val = supremum(I);
        ksi = sign(temp.G)';
    elseif strcmp(type,'lower')
        val = infimum(I); 
        ksi = -sign(temp.G)';
    elseif strcmp(type,'range')
        val = I;
        ksi = [-sign(temp.G)', sign(temp.G)'];
    end

    % calculate support vector
    if nargout >= 2
        x = cZ.c + cZ.G*ksi;
    end
    return

end

% object properties
A = cZ.A;
b = cZ.b;
n = size(A,2);

% ksi in [-1, 1]
lb = -ones(n,1);
ub = ones(n,1);

% project the zonotope onto the direction
f = dir' * cZ.G;

% init values for better code logic
fval = [];
x = [];
ksi = [];

% special cases
if ~any(f)
    % projection of the zonotope equals 0 -> no extension of the set in
    % this direction
    fval = 0;
    ksi = zeros(n,1);

elseif nargout == 1

    % check if the objective function is identical to the normal vector of
    % a constraint -> only one solution
    idxPos = all(withinTol(f,A),2);
    idxNeg = all(withinTol(-f,A),2);
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

persistent isMosek
if isempty(isMosek)
    isMosek = isSolverInstalled('mosek');
end

if isempty(fval)
    % linear program:
    % max_{x \in cZ} dir' * x
    % s.t. x = c + G*ksi
    %      A * ksi = b
    
    if isMosek

        % lower bound
        if any(strcmp(type,{'lower','range'}))
            res = msklpopt(f',A,b,b,lb,ub,[],'minimize echo(0)');

            if strcmp(res.sol.itr.prosta,'PRIMAL_INFEASIBLE')
                % primal infeasible -> empty set
                if strcmp(type,'range')
                    val = interval.empty(1);
                elseif strcmp(type,'lower')
                    val = Inf;
                end
                return
            elseif strcmp(res.sol.itr.prosta,'DUAL_INFEASIBLE')
                % should not be unbounded -> solver issue
                throw(CORAerror('CORA:solverIssue'));
            elseif strcmp(res.sol.itr.prosta,'PRIMAL_AND_DUAL_FEASIBLE')
                % combine factors
                ksi = res.sol.itr.xx;
                fval = res.sol.itr.dobjval;
            end

        end

        % upper bound
        if any(strcmp(type,{'upper','range'}))
            res = msklpopt(f',A,b,b,lb,ub,[],'maximize echo(0)');

            if strcmp(res.sol.itr.prosta,'PRIMAL_INFEASIBLE')
                % primal infeasible -> empty set
                if strcmp(type,'range')
                    val = interval.empty(1);
                elseif strcmp(type,'upper')
                    val = -Inf;
                end
                return
            elseif strcmp(res.sol.itr.prosta,'DUAL_INFEASIBLE')
                % should not be unbounded -> solver issue
                throw(CORAerror('CORA:solverIssue'));
            elseif strcmp(res.sol.itr.prosta,'PRIMAL_AND_DUAL_FEASIBLE')
                % combine factors
                ksi = [ksi, res.sol.itr.xx];
                fval = [fval, res.sol.itr.dobjval];
            end
        end

    else
        % MATLAB linprog

        % linear program options
        persistent options
        if isempty(options)
            options = optimoptions('linprog','Algorithm','dual-simplex','display','off');
        end

        % lower bound
        if any(strcmp(type,{'lower','range'}))
            [ksi,fval,exitflag] = linprog(f',[],[],A,b,lb,ub,options);

            if exitflag == -2
                % primal infeasible -> empty set
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
            [ksi_,fval_,exitflag] = linprog(-f',[],[],A,b,lb,ub,options);

            if exitflag == -2
                % primal infeasible -> empty set
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
end

% calculate bound by adding the zonotope center
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

% ------------------------------ END OF CODE ------------------------------
