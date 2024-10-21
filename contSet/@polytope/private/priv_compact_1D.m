function [A,b,Ae,be,empty] = priv_compact_1D(A,b,Ae,be,tol)
% priv_compact_1D - removes all redundant constraints for a 1D polytope
%
% Syntax:
%    [A,b,Ae,be,empty] = priv_compact_1D(A,b,Ae,be,tol)
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

if ~isempty(Ae)
    % normalization yields equality constraints of the form
    %   ax = 0, ax = 1
    % normalize with respect to b
    [A,b,Ae,be] = priv_normalizeConstraints(A,b,Ae,be,'b');

    % more than one unique equality constraint is not feasible
    if size(Ae,1) > 1
        if any(withinTol(be,1,tol)) && any(withinTol(be,0,tol))
            % any combination Ax = 1 and Ax = 0 -> empty
            empty = true; A = []; b = []; Ae = []; be = [];
        end
        
        % indices where b = 1 and b = 0
        idx_1 = withinTol(be,1);
        Ae_1 = Ae(idx_1);
        idx_0 = withinTol(be,0);
        Ae_0 = Ae(idx_0);

        if nnz(idx_1) > 1 && ~all(withinTol(Ae_1,Ae_1(1),tol))
            % more than one constraint with ... = 1
            % and not the same value in A -> empty
            empty = true; A = []; b = []; Ae = []; be = [];
            return
        elseif nnz(idx_0) > 1 && ~all(withinTol(Ae_0,Ae_0(1),tol))
            % more than one constraint with ... = 0
            % and not the same value in A -> empty
            empty = true; A = []; b = []; Ae = []; be = [];
            return
        end

        % all equality constraints are the same as the first one; take
        % this one and normalize w.r.t Ae
        Ae = 1; be = be(1) / Ae(1);
    end
end

if ~isempty(A)
    % normalize rows of A matrix
    [A,b,~,be_] = priv_normalizeConstraints(A,b,Ae,be,'A');

    % take outermost halfspaces
    idxPos = A > 0;
    A = [1;-1];
    b = [min(b(idxPos));min(b(~idxPos))];

    % check if empty (including equality constraint)
    if length(b) < 2
        % only one inequality constraint
        if all(idxPos)
            A = 1;
        elseif ~all(idxPos)
            A = -1;
        end
    elseif b(1) < -b(2)
        % constraints of the form
        %   ax <= b, ax >= b_, where b_ > b
        % -> infeasible!
        empty = true; A = []; b = []; Ae = []; be = [];
        return
    elseif ~isempty(be)
        % there are equality constraints
        if any(be_ > b(1)) || any(be_ < -b(2))
            % additionally, an equality constraint that does not comply
            % with the inequality constraints
            empty = true; A = []; b = []; Ae = []; be = [];
            return
        else
            % equality constraint is satisfied by inequality
            % constraints, only keep equality constraint
            A = []; b = [];
        end
    end
    
end

% check if inequalities and equality are feasible
% note: at most 2 inequalities and 1 equality
% if ~isempty(A) && ~isempty(Ae)
%     % all are normalized w.r.t to A, so check b values
%     if A(1) == 1 && b(1) 
% 
% end

% ------------------------------ END OF CODE ------------------------------
