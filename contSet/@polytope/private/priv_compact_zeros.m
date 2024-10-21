function [A,b,Ae,be,empty] = priv_compact_zeros(A,b,Ae,be,tol)
% priv_compact_zeros - removes all constraints where the vector is all-zero
%    1) removes all all-zero inequality constraints, i.e., 0*x <= a
%    also, if a < 0, then the constraint is infeasible and the set is empty
%    2) removes all all-zero equality constraints, i.e., 0*x = a, a = 0
%    also, if a ~= 0, then the constraint is infeasible and the set is empty
%
% Syntax:
%    [A,b,Ae,be,empty] = priv_compact_zeros(A,b,Ae,be,tol)
%
% Inputs:
%    A - inequality constraint matrix
%    b - inequality constraint offset
%    Ae - equality constraint matrix
%    be - equality constraint offset
%    tol - tolerance
%
% Outputs:
%    A - inequality constraint matrix
%    b - inequality constraint offset
%    Ae - equality constraint matrix
%    be - equality constraint offset
%    empty - true/false whether polytope is empty
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

empty = false;

% remove trivially redundant rows with 0*x <= b
A_zero_rows = find(all(withinTol(A,0),2));
if any(A_zero_rows)
    % check whether any of those constraints is infeasible
    if any(b(A_zero_rows) < 0 & ~withinTol(b(A_zero_rows),0,tol))
        % constraint 0*x <= b with b < 0  =>  infeasible, i.e., empty set
        empty = true; A = []; b = []; Ae = []; be = [];
        return
    else
        % all constraints are 0*x <= b where all b_i > 0
        A(A_zero_rows,:) = [];
        b(A_zero_rows) = [];
    end
end

% remove trivially redundant rows with 0*x = 0
Ae_zero_rows = find(all(withinTol(Ae,0),2));
if any(Ae_zero_rows)
    % check whether any of those constraints is infeasible
    if any(~withinTol(be(Ae_zero_rows),0,tol))
        % constraint 0*x = be with be ~= 0  =>  infeasible, i.e., empty set
        empty = true; A = []; b = []; Ae = []; be = [];
        return
    else
        % all constraints are 0*x = 0
        Ae(Ae_zero_rows,:) = [];
        be(Ae_zero_rows) = [];
    end
end

% ------------------------------ END OF CODE ------------------------------
