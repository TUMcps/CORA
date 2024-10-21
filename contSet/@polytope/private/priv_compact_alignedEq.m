function [Ae,be,empty] = priv_compact_alignedEq(Ae,be,tol)
% priv_compact_alignedEq - removes all redundant aligned constraints
%    note: expects normalized constraints with respect to 'A'
%
% Syntax:
%    [Ae,be,empty] = priv_compact_alignedEq(Ae,be,tol)
%
% Inputs:
%    Ae - equality constraint matrix
%    be - equality constraint offset
%    tol - tolerance
%
% Outputs:
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

% assume resulting set is non-empty
empty = false;

% number of equality constraints and dimension
[nrConEq,n] = size(Ae);

% try to find aligned constraint
% vectors -> either infeasible (different value for be) or
% redundant (same value for be)

% pre-compute dot product
dotprod_norm = Ae * Ae';

% init logical array for which constraints to check
irredundantEq = true(size(Ae,1),1);

for i=1:nrConEq
    % if there are aligned constraints, then they have to have the
    % same value in be (all but one of them are redundant),
    % otherwise the polytope is the empty set because there is no
    % point that can satisfy to parallel Ae*x = be at the same time

    % only check if the i-th equality constraint has not already
    % been removed due to alignment
    if irredundantEq(i)
        % check for aligned vectors
        alignedConstraints = all(...
            withinTol(Ae(i,:) - dotprod_norm(:,i)*Ae(i,:),...
            zeros(1,n)), 2);
    
        if nnz(alignedConstraints) > 1
            % at least two constraints are aligned
            if all(withinTol(be(alignedConstraints),be(i)))
                % remove all constraints but the first one
                irredundantEq(alignedConstraints) = false;
                irredundantEq(i) = true;
            else
                % polytope is the empty set
                empty = true; return
            end
        end
    end
end

% remove all redundant equality constraints
Ae = Ae(irredundantEq,:);
be = be(irredundantEq);

% ------------------------------ END OF CODE ------------------------------
