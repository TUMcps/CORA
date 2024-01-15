function [val,x] = supportFunc_(P,dir,type,varargin)
% supportFunc_ - Calculate the upper or lower bound of a polytope object
%    along a certain direction
%
% Syntax:
%    val = supportFunc_(P,dir,type)
%    [val,x] = supportFunc_(P,dir,type)
%
% Inputs:
%    P - polytope object
%    dir - direction for which the bounds are calculated (vector of size
%          (n,1) )
%    type - upper or lower bound ('lower', 'upper', 'range')
%
% Outputs:
%    val - bound of the constraind zonotope in the specified direction
%    x - support vector
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/supportFunc, conZonotope/supportFunc_

% Authors:       Niklas Kochdumper, Victor Gassmann, Mark Wetzlinger
% Written:       19-November-2019
% Last update:   16-March-2021 (added unbounded support)
%                10-December-2022 (MW, add 'range')
%                13-December-2022 (MW, add call to MOSEK)
%                15-November-2023 (MW, computation for vertex representation)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check fullspace
if representsa_(P,'fullspace',0)
    x = Inf .* sign(dir);
    x(isnan(x)) = 0;
    switch type
        case 'upper'
            val = Inf;
        case 'lower'
            val = -Inf;
        case 'range'
            val = interval(-Inf,Inf);
    end
    return
end

% check if vertex representation given (skips linear program)
if ~isempty(P.V.val)
    if strcmp(type,'upper')
        vals = dir' * P.V.val;
        [val,idx] = max(vals);
        x = P.V.val(:,idx);
    elseif strcmp(type,'lower')
        vals = -dir' * P.V.val;
        [val,idx] = min(-vals);
        x = P.V.val(:,idx);
    elseif strcmp(type,'range')
        vals_upper = dir' * P.V.val;
        vals_lower = -dir' * P.V.val;
        [val_upper,idx_upper] = max(vals_upper);
        [val_lower,idx_lower] = min(-vals_lower);
        val = interval(val_lower, val_upper);
        x = [P.V.val(:,idx_lower), P.V.val(:,idx_upper)];
    end
    return
end

% check if MOSEK is installed
persistent isMosek
if isempty(isMosek)
    isMosek = isSolverInstalled('mosek');
end

% linear program options (only if MATLAB linprog is used)
persistent options
if isempty(options)
    options = optimoptions('linprog','display','off');
end

% upper or lower bound
if strcmp(type,'lower') || strcmp(type,'upper')
    [val,x] = aux_solveLinProg(P,dir,type,isMosek,options);
elseif strcmp(type,'range')
    [val_upper,x_upper] = aux_solveLinProg(P,dir,'upper',isMosek,options);
    [val_lower,x_lower] = aux_solveLinProg(P,dir,'lower',isMosek,options);
    % combine values
    val = interval(val_lower,val_upper);
    x = [x_lower x_upper];
end

end


% Auxiliary functions -----------------------------------------------------

function [val,x] = aux_solveLinProg(P,dir,type,isMosek,options)

if strcmp(type,'upper')
    s = -1;
elseif strcmp(type,'lower')
    s = 1;
end

% simple check: empty polytope
if isempty(P.A) && isempty(P.Ae)
    val = s*Inf; x = [];
    return
end

% solve linear program
if isMosek

    % number of constraints
    nrCon = length(P.b);

    % rewrite for MOSEK syntax
    c = dir';
    a = [P.A; P.Ae];
    blc = [-Inf(nrCon,1); P.be];
    buc = [P.b; P.be];
    
    % call MOSEK
    if strcmp(type,'upper')
        res = msklpopt(c,a,blc,buc,[],[],[],'maximize echo(0)');
    elseif strcmp(type,'lower')
        res = msklpopt(c,a,blc,buc,[],[],[],'minimize echo(0)');
    end

    % read out value of support function and support vector
    if strcmp(res.sol.itr.prosta,'PRIMAL_AND_DUAL_FEASIBLE')
        val = res.sol.itr.dobjval;
        x = res.sol.itr.xx;
    else
        % either unbounded or infeasible
        if strcmp(res.sol.itr.prosta,'DUAL_INFEASIBLE')
            % unbounded
            % TODO: check value for x (does not have to be Inf everywhere!)
            val = -s*Inf;
            x = -s*sign(dir).*Inf(length(dir),1);
        elseif strcmp(res.sol.itr.prosta,'PRIMAL_INFEASIBLE')
            % infeasible -> empty set
            val = s*Inf;
            x = [];
        end
    end

else
    
    % solve linear program
    [x,val,exitflag] = linprog(s*dir',P.A,P.b,P.Ae,P.be,[],[],options);
    val = s*val;

    if exitflag == -3
        % unbounded
        val = -s*Inf;
        x = -s*sign(dir).*Inf(length(dir),1);
    elseif exitflag == -2
        % infeasible -> empty set
        val = s*Inf;
        x = [];
    elseif exitflag ~= 1
        throw(CORAerror('CORA:solverIssue'));
    end
end

end

% ------------------------------ END OF CODE ------------------------------
