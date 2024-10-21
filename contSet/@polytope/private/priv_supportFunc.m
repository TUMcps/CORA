function [val,x] = priv_supportFunc(A,b,Ae,be,dir,type)
% priv_supportFunc - computes the halfspace representation of the box enclosure
%    given a vertex representation
%
% Syntax:
%    [val,x] = priv_supportFunc(A,b,Ae,be,dir,type)
%
% Inputs:
%    A - inequality constraint matrix
%    b - inequality constraint offset
%    Ae - equality constraint matrix
%    be - equality constraint offset
%    dir - direction
%    type - 'upper' or 'lower'
%
% Outputs:
%    val - value of the support function
%    x - support vector
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       03-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if strcmp(type,'upper')
    s = -1;
elseif strcmp(type,'lower')
    s = 1;
end

% simple check: empty polytope (fullspace)
if isempty(A) && isempty(Ae)
    val = -s*Inf; x = [];
    return
end

% set up linear program
problem.f = s*dir';
problem.Aineq = A;
problem.bineq = b;
problem.Aeq = Ae;
problem.beq = be;
problem.lb = [];
problem.ub = [];

% solve linear program
[x,val,exitflag] = CORAlinprog(problem);
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

% ------------------------------ END OF CODE ------------------------------
