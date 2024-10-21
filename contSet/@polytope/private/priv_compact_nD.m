function [A,b,empty] = priv_compact_nD(A,b,Ae,be,n,tol)
% priv_compact_nD - removes all redundant inequality constraints in the
%    representation of an nD polytope
%
% Syntax:
%    [A,b,empty] = priv_compact_nD(A,b,Ae,be,n,tol)
%
% Inputs:
%    A - inequality constraint matrix
%    b - inequality constraint offset
%    Ae - equality constraint matrix
%    be - equality constraint offset
%    n - dimension of the polytope
%    tol - tolerance
%
% Outputs:
%    A - inequality constraint matrix
%    b - inequality constraint offset
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
if isempty(b)
    return
end

% remove redundant inequality constraints
[A,b,irredundantIneq,done] = aux_irredundantIneq(A,b);
if done
    return
end

% check other heuristics, e.g.,
% - any inequality that contains all the extreme points of the bounding
%   box is redundant (caution: bounding box also uses linprog -> could
%   only be useful if number of halfspaces >> dim(P))

% all code until here attempts to avoid solving linear programs as much
% as possible; we now perform a single linear program upfront to check
% for emptiness
if priv_representsa_emptySet(A,b,Ae,be,n,eps)
    empty = true; A = []; b = [];
    return
end

% for remaining rows in inequality constraint matrix: test if i-th
% inequality is non-redundant by computing the support function for the
% polytope without the i-th inequality in the direction of the i-th
% inequality and comparing the result to the corresponding offset

% update number of constraints
nrConIneq = size(A,1);

% init indices for constraints that should be kept
idxKeep = true(nrConIneq,1);
% index for temporarily-kept constraints
idxCurrKept = true(nrConIneq,1);

% loop over all remaining inequality constraints
for i=1:nrConIneq
    % solve primal LP for polytope without i-th constraint

    % skip computation if constraint has already been determined to be
    % irredundant by heuristics above or in the computation below
    if ~irredundantIneq(i)
    
        % remove i-th constraint
        idxCurrKept(i) = false;
        H = A(idxKeep & idxCurrKept,:);
        d = b(idxKeep & idxCurrKept);

        % compute support function in the direction of the i-th constraint
        [val,extreme_point] = priv_supportFunc(H,d,Ae,be,A(i,:)','upper');
        
        % reduced polytope is empty -> original polytope is, too
        if val == -Inf
            empty = true; A = []; b = [];
            return
        end

        if val < b(i) || withinTol(val,b(i),tol)
            % redundant
            idxKeep(i) = false;
            % check which constraint are active for that extreme point
            % -> they cannot be redundant
            cons_active = withinTol(A*extreme_point,b,tol);
            irredundantIneq = irredundantIneq | cons_active;
        end
    
        % go to next constraint (does not matter if this constraint has
        % just been removed, because idxKeep will not allow it to be
        % used anymore)
        idxCurrKept(i) = true;

    end
end

% remove redundant inequality constraints
A = A(idxKeep,:);
b = b(idxKeep);

end


% Auxiliary functions -----------------------------------------------------

function [A,b,irredundantIneq,done] = aux_irredundantIneq(A,b)

% number of inequality constraints and dimension
[nrConIneq,n] = size(A);

% heuristics: find orthogonal and aligned vectors
irredundantIneq = false(size(A,1),1);
redundantIneq = false(size(A,1),1);

% pre-compute dot product
dotprod_norm = A * A';

for i=1:nrConIneq
    % a constraint is orthogonal if dot product with all other
    % constraints is 0 or -1 -> then this constraint is irredundant;
    % conversely, if there are aligned constraints, then the one with
    % the larger value for b in Ax <= b is definitely redundant

    % check for aligned vectors
    alignedConstraints = all(...
        withinTol(A(i,:) - dotprod_norm(:,i)*A(i,:),...
        zeros(1,n)), 2);

    if nnz(alignedConstraints) > 1
        % check which one has the smallest offset (b)
        bmin = min(b(alignedConstraints));
        % logical array for all aligned vectors with larger b -> add to
        % redundancies
        redundantIneq = redundantIneq ...
            | xor(b == bmin & alignedConstraints,alignedConstraints);
    end

    % check if all other constraints are orthogonal to this one
    if all(withinTol(dotprod_norm(~alignedConstraints,i),0))
        irredundantIneq(i) = true;
    end
end
% remove redundant constraints following alignment criterion
A(redundantIneq,:) = [];
b(redundantIneq) = [];
% adapt indices of list of irredundant constraints
irredundantIneq = irredundantIneq(~redundantIneq);

% check if minimal representation already obtained
done = all(irredundantIneq);

end

% ------------------------------ END OF CODE ------------------------------
