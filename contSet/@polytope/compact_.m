function P_out = compact_(P,type,tol,varargin)
% compact_ - Computes the irredundant H/V-representation of a polytope;
%    inequality constraints: the main idea is to compare the offset of a
%    given constraint with the support function of the polytope without
%    that constraints; specialized methods for 1D and 2D polytopes
%
% Syntax:
%    P_out = compact_(P)
%    P_out = compact_(P,type)
%    P_out = compact_(P,type,tol)
%
% Inputs:
%    P - polytope object
%    type - (optional) set of constraints that should be reduced
%           'all' (default): inequality and equality constraints
%           'A': only inequality constraints
%           'Ae': only equality constraints
%           'aligned': check only for constraints with same orientation
%           'zeros': check for 0*x <= a and 0*x = a constraints
%           'V': remove redundancies in vertex representation
%
% Outputs:
%    P_out - polytope object in reduced or minimal representation
%
% Example:
%    A = [1 1; 1 0; 1 1; -1 0; 0 -1; 1 1; 1 2];
%    b = [7; 3; 5; -1; -1; 6; 20];
%    P = polytope(A,b);
%    P_ = compact(P);
%
%    figure; hold on;
%    plot(P);
%    plot(P_,[1,2],'r--');
%
% Reference: MPT-Toolbox https://www.mpt3.org/
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/compact

% Authors:       Viktor Kotsev, Mark Wetzlinger
% Written:       23-May-2022
% Last update:   08-December-2022 (MW, use logical indexing to eliminate inequalities in sequential tests)
%                09-December-2022 (MW, pre-processing using orthogonality and alignment criteria)
%                12-December-2022 (MW, fix faulty alignment criterion)
%                14-December-2022 (MW, add support for MOSEK)
%                04-April-2023 (MW, use persistent variables)
%                15-January-2024 (TL, bug fix single infeasible equality constraint)
% Last revision: 28-May-2023 (MW, increase readability via aux_ functions, integrate code from removeRedundancies)
%                31-July-2023 (MW, rename 'compact_', integrate deleteZeros, integrate minVRep)

% ------------------------------ BEGIN CODE -------------------------------

% quickest exit
if representsa_(P,'fullspace',0)
    P_out = polytope.Inf(dim(P));
    return
end

% copy set (also copies properties)
P_out = polytope(P);

% reduce vertex representation
if strcmp(type,'V') || (isempty(P.A) && isempty(P.Ae) && ~isempty(P.V.val))
    P_out = aux_minVRep(P_out,tol);
    % the set with the reduced representation size has the same properties
    % (boundedness, etc.) as the original one
    P = aux_copyProperties(P,P_out);
    return
end

% remove all-zero constraints
P_out = aux_deleteZeros(P_out,tol);
% exit if this is all that was required
if strcmp(type,'zeros')
    P = aux_copyProperties(P,P_out);
    return
end

% quick exits:
if length(P_out.b) == 1 && isempty(P_out.be)
    % only one inequality, no equality given -> already minimal
    P_out = aux_oneHalfspace(P_out,tol);
    P = aux_copyProperties(P,P_out);
    return
elseif isempty(P_out.b) && length(P_out.be) == 1
    % only one equality, no inequality given -> already minimal
    P_out.minHRep.val = true;
    P = aux_copyProperties(P,P_out);
    return
end

% normalize rows for all subsequent operations
P_out = normalizeConstraints(P_out,'A');

% removal of aligned constraints
P_out = aux_removeAlignedIneq(P_out);
if strcmp(type,'aligned') || (length(P_out.b) <= 1 && isempty(P_out.be))
    P = aux_copyProperties(P,P_out);
    return
end

% quick check for two inequality constraints
if length(P_out.b) == 2 && isempty(P_out.Ae)
    P_out = aux_twoHalfspaces(P_out);
    P_out.minHRep.val = true;
    % check whether object has become fully empty
    P = aux_copyProperties(P,P_out);
    return
end

% dimension
n = dim(P_out);

% special algorithms for 1D and 2D
if n == 1
    P_out = aux_1D(P_out,type); 
    P = aux_copyProperties(P,P_out);
    return
elseif n == 2
    P_out = aux_2D(P_out,type); 
    P_out.minHRep.val = true;
    P = aux_copyProperties(P,P_out);
    return
end

% read out properties
A = P_out.A; b = P_out.b;
Ae = P_out.Ae; be = P_out.be;

% remove redundant equality constraints
if ~isempty(Ae) && (strcmp(type,'all') || strcmp(type,'Ae'))

    [Ae,be,empty] = aux_removeEqCon(Ae,be);

    % quick exit if two equality constraints cannot be fulfilled at the
    % same time
    if empty
        P_out = polytope.empty(n);
        P = aux_copyProperties(P,P_out);
        return
    end
end


% remove redundant inequality constraints
if ~isempty(A) && (strcmp(type,'all') || strcmp(type,'A'))

    [A,b,irredundantIneq,done] = aux_irredundantIneq(A,b);
    if done
        P_out.A = A; P_out.b = b;
        P_out.minHRep.val = true;
        P = aux_copyProperties(P,P_out);
        return
    end
    
    % check other heuristics, e.g.,
    % - any inequality that contains all the extreme points of the bounding
    %   box is redundant (caution: bounding box also uses linprog -> could
    %   only be useful if number of halfspaces >> dim(P))
    
    % all code until here attempts to avoid solving linear programs as much
    % as possible; we now perform a single linear program upfront to check
    % for emptiness
    if representsa_(polytope(A,b,Ae,be),'emptySet',eps)
        P_out = polytope.empty(n);
        P = aux_copyProperties(P,P_out);
        return
    end
    
    % linear programming
    persistent isMosek
    if isempty(isMosek)
        isMosek = isSolverInstalled('mosek');
    end
    persistent options
    if isempty(options)
        options = optimoptions('linprog','display','off');
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
        % irredundant by heuristics above
        if ~irredundantIneq(i)
        
            % remove i-th constraint
            idxCurrKept(i) = false;
            H = A(idxKeep & idxCurrKept,:);
            d = b(idxKeep & idxCurrKept);
            
            if isMosek
                
                % rewrite to MOSEK syntax
                a = [H; Ae];
                blc = [-Inf(size(H,1),1); be];
                buc = [d; be];
    
                % call MOSEK
                res = msklpopt(A(i,:)',a,blc,buc,[],[],[],'maximize echo(0)');
    
                if strcmp(res.sol.itr.prosta,'PRIMAL_AND_DUAL_FEASIBLE')
                    % value of support function in that direction
                    minvalue = A(i,:) * res.sol.itr.xx;
                    if minvalue < b(i) || withinTol(minvalue,b(i),tol)
                        % redundant
                        idxKeep(i) = false;
                    end
                elseif strcmp(res.sol.itr.prosta,'PRIMAL_INFEASIBLE')
                    % reduced set of constraints is infeasible -> empty set
                    P_out = polytope.empty(n);
                    P = aux_copyProperties(P,P_out);
                    return
                end
    
            else
                % MATLAB linprog
    
                % solve LP
                [x,~,exitflag] = linprog(-A(i,:),H,d,Ae,be,[],[],options);
            
                % solution has been found
                if exitflag > 0
                    minvalue = A(i,:)*x;
                    if minvalue < b(i) || withinTol(minvalue,b(i),tol)
                        % redundant
                        idxKeep(i) = false;
                    end
                elseif exitflag == -2
                    % reduced set of constraints is infeasible -> empty set
                    P_out.A = zeros(0,n); P_out.b = [];
                    P_out.Ae = zeros(0,n); P_out.be = [];
                    P_out.minHRep.val = true;
                    P_out.emptySet.val = true;
                    P_out.bounded.val = true;
                    P_out.fullDim.val = false;
                    P_out.V.val = zeros(n,0);
                    P_out.minVRep.val = true;
                    P = aux_copyProperties(P,P_out);
                    return
                elseif exitflag == -4
                    throw(CORAerror('CORA:solverIssue'));
                end
            end
        
            % go to next constraint (does not matter if this constraint has
            % just been removed, because idxKeep will not allow it to be
            % used anymore)
            idxCurrKept(i) = true;
    
        end
    end
    
    % remove redundant inequality constraints
    P_out = polytope(A(idxKeep,:), b(idxKeep), Ae, be);
    P_out.minHRep.val = true;

end

end


% Auxiliary functions -----------------------------------------------------

function P = aux_copyProperties(P,P_out)
% copies all properties that were determined during the redundancy removal
% from P_out back to P
% note: at the beginning, all properties were copied from P to P_out, so
% assigning them back does not override anything

    % full-dimensionality, bounded, emptiness do not depend on any redudancy
    % in the representation
    P.emptySet.val = P_out.emptySet.val;
    P.fullDim.val = P_out.fullDim.val;
    P.bounded.val = P_out.bounded.val;
    
    % P is only in minimal representation, if it has the same number of
    % inequality and equality constraints as P_out
    P.minHRep.val = length(P_out.b) == length(P.b) ...
        && length(P_out.be) == length(P.be);
    
    % copy vertices as well
    P.V.val = P_out.V.val;
    P.minVRep.val = P_out.minVRep.val;

end

function P = aux_minVRep(P,tol)
% removes redundancies in the vertex representation

    % check for V-representation
    if isempty(P.V.val)
        % nothing to do here...
        return;
    end
    
    % extract arguments
    V = P.V.val;
    
    % handle 1-dimensional case separately
    if dim(P) == 1
        min_V = min(min(V));
        max_V = max(max(V));
        V = [min_V; max_V]';    
        P.V.val = V;
    else
        % compute convex hull and take only the unique points
        try
            K = convhulln(V');
            P.V.val = V(:,unique(K));
        catch ME
            if strcmp(ME.identifier,'MATLAB:cgprechecks:NotEnoughPts')
                % not enough unique points specified -> we assume that this
                % means minimal representation
            else
                rethrow(ME);
            end
        end
    end
    
    % set minVRep to true
    P.minVRep.val = true;

end

function P = aux_deleteZeros(P,tol)
% 1) removes all all-zero inequality constraints, i.e., 0*x <= a
% also, if a < 0, then the constraint is infeasible and the set is empty
% 2) removes all all-zero equality constraints, i.e., 0*x = a, a = 0
% also, if a ~= 0, then the constraint is infeasible and the set is empty
% note: this function is always called! (except reduction of V-rep.)

    % dimension
    n = dim(P);

    % remove trivially redundant rows with 0*x <= b
    A_zero_rows = find(all(withinTol(P.A,0),2));
    if any(A_zero_rows)
        % check whether any of those constraints is infeasible
        if any(P.b(A_zero_rows) < 0 & ~withinTol(P.b(A_zero_rows),0,tol))
            % constraint 0*x <= b with b < 0  =>  infeasible, i.e., empty set
            P = polytope.empty(n);
            return
        else
            % all constraints are 0*x <= b where all b_i > 0
            P.A(A_zero_rows,:) = [];
            P.b(A_zero_rows) = [];
        end
    end

    % remove trivially redundant rows with 0*x = 0
    Ae_zero_rows = find(all(withinTol(P.Ae,0),2));
    if any(Ae_zero_rows)
        % check whether any of those constraints is infeasible
        if any(~withinTol(P.be(Ae_zero_rows),0,tol))
            % constraint 0*x = be with be ~= 0  =>  infeasible, i.e., empty set
            P = polytope.empty(n);
            return
        else
            % all constraints are 0*x = 0
            P.Ae(Ae_zero_rows,:) = [];
            P.be(Ae_zero_rows) = [];
        end
    end

    % if no constraints are left, but none have been infeasible, we have an
    % unbounded set (fullspace)
    if isempty(P.A) && isempty(P.Ae)
        P.emptySet.val = false;
        P.bounded.val = false;
        P.fullDim.val = true;
    end    

end

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

function [Ae,be,empty] = aux_removeEqCon(Ae,be)

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

end

function P = aux_2D(P,type)
% special computation for 2D case without linear programs:
% First, we order the (already normalized) constraints by angle. Then, we
% check whether the middle constraint between two adjacent constraints is
% required by computing the intersecting vertex between the outer two
% constraints and checking whether that middle constraint has a smaller
% support function value than the intersecting vertex in the direction
% of the middle constraint

    n = dim(P);
    Ae = P.Ae;
    be = P.be;

    if ~isempty(P.Ae) && (strcmp(type,'all') || strcmp(type,'Ae'))
        % remove duplicates
        [Ae,be,empty] = aux_removeEqCon(P.Ae,P.be);

        if empty
            % aux_removeEqCon already determined emptiness
            P = polytope.empty(2);
            return
        end

    end

    if ~isempty(P.A) && (strcmp(type,'all') || strcmp(type,'A'))

        % transpose constraints, read out number
        A = P.A';
        nrCon = length(P.b);
    
        % constraints are normalized, order by angle
        angles = atan2d(A(2,:),A(1,:));
        [angles,idx] = sort(angles);
        A = A(:,idx);
        b = P.b(idx);
    
        % remove parallel constraints (need to be in order!)
        idxKept = true(1,nrCon);
        startIdx = nrCon; endIdx = 1; idxUnchecked = true(1,nrCon);
        
        % loop until all have been compared once
        while any(idxUnchecked)
    
            % compute dot product, check for alignment
            if idxUnchecked(startIdx) && withinTol(A(:,startIdx)' * A(:,endIdx),1)
                % compute dot product with all constraints, compare all which
                % are parallel and keep the one with the smallest offset value
                dotprod = A(:,startIdx)' * A;
                idxParallel = withinTol(dotprod,1);
                [~,minIdx] = min(b(idxParallel));
    
                % remove all with larger offset
                idxKept = idxKept & ~idxParallel;
    
                % put back the one with smallest offset value
                idxKeep = find(idxParallel,minIdx);
                idxKept(idxKeep(end)) = true;
    
                % all from this batch have been checked -> do not check again
                idxUnchecked(idxParallel) = false;
            end
    
            % increment indices
            idxUnchecked(startIdx) = false;
            startIdx = endIdx;
            endIdx = endIdx + 1;
        end
        A = A(:,idxKept);
        b = b(idxKept);
        angles = angles(idxKept);
        nrCon = length(b);
    
        % if only two halfspaces are left
        if length(b) == 2
            P = aux_twoHalfspaces(polytope(A',b));
            P.Ae = Ae; P.be = be;
            return
        end
    
        % 90 degree counter-clockwise shifted halfspaces
        A_ = [A(2,:); -A(1,:)];
    
        % indices for constraints that are kept
        idxKept = true(1,nrCon);
        % indices for constraints hhat have been checked
        idxUnchecked = true(1,nrCon);
    
        startIdx = 1; middleIdx = 2; endIdx = 3;
    %     make_plot = false;
    %     if make_plot
    %     figure;
    %     end
        while any(idxUnchecked)
    
            % halfspace at middleIdx will now be checked
            idxUnchecked(middleIdx) = false;
    
            % for potential redundancy, angle between startIdx and endIdx must
            % not be larger than 180 degrees, otherwise intersection at the
            % wrong side
            if aux_anglesTooLarge(angles(startIdx),angles(endIdx))
                % increment indices (take care of wrap-around!)
                [startIdx,middleIdx,endIdx] = ...
                    aux_incrementIndices(startIdx,middleIdx,endIdx,nrCon);
                continue
            end
    
            % check if constraint between startIdx and endIdx is necessary by
            % computing the intersection of
            %   p1 + a*p2vec, a \in R
            %   p2 + b*p2vec, b \in R
            % which is
            %   p1 + a*p2vec = p2 + b*p2vec
            %   a*p1vec - b*p2vec = p2 - p1
            %   [p1vec -p2vec] * [a;b] = p2 - p1
            p1 = A(:,startIdx) * b(startIdx);
            p2 = A(:,endIdx) * b(endIdx);
            p1vec = A_(:,startIdx);
            p2vec = -A_(:,endIdx);
            
            % if vectors are parallel, then middleIdx non-redundant
            if rank([p1vec p2vec]) < 2
                [startIdx,middleIdx,endIdx] = ...
                    aux_incrementIndices(startIdx,middleIdx,endIdx,nrCon);
                continue;
            end
            
            factors = [p1vec p2vec] \ (p2-p1);
            % plug in
            x = p1 + factors(1)*A_(:,startIdx);
            % ...should be the same as p2+factors(2)*(-A_(:,endIdx))
    
            % debug ---
    %         if make_plot
    %         scatter(0,0,16,'r','filled');
    %         hold on; box on; axis equal
    %         plot(interval([-1;-1],[1;1]),[1,2],'r');
    %         scatter(p1(1),p1(2),16,'k','filled');
    %         scatter(p2(1),p2(2),16,'k','filled');
    %         line1 = [p1-10*p1vec p1+10*p1vec];
    %         line2 = [p2-10*p2vec p2+10*p2vec];
    %         plot(line1(1,:),line1(2,:),'b');
    %         plot(line2(1,:),line2(2,:),'b');
    %         scatter(x(1),x(2),16,'green','filled');
    %         p3 = A(:,middleIdx) * b(middleIdx);
    %         p3vec = A_(:,middleIdx);
    %         line3 = [p3-10*p3vec p3+10*p3vec];
    %         plot(line3(1,:),line3(2,:),'r');
    %         hold off;
    %         end
            % debug ---
    
            % evaluate support function of intersecting point and compare to
            % offset of middleIdx'th halfspace
            val = A(:,middleIdx)' * x;
            if b(middleIdx) > val || withinTol(b(middleIdx),val,1e-9)
                % -> middle halfspace is redundant
                idxKept(middleIdx) = false;
                % increment indices (take care of wrap-around!)
                middleIdx = endIdx;
                endIdx = endIdx+1;
                if endIdx > nrCon; endIdx = 1; end
            else
                % -> middle halfspace is non-redundant
                % increment indices (take care of wrap-around!)
                [startIdx,middleIdx,endIdx] = ...
                    aux_incrementIndices(startIdx,middleIdx,endIdx,nrCon);
            end
    
        end
    
        % remove halfspaces
        P.A = A(:,idxKept)';
        P.b = b(idxKept);
        P.Ae = Ae;
        P.be = be;

        % TODO: check if empty without linprog (?)
    end

end

function res = aux_anglesTooLarge(alpha,beta)
% auxiliary function for 2D case
% checks if angle between alpha and beta is larger than 180 degrees
% alpha, beta \in [-180,180]

    % shift by 180
    alpha = alpha + 180; beta = beta + 180;
    
    % shift alpha and beta such that alpha is at 0 (alpha not needed)
    beta = mod(beta + (360-alpha),360);
    
    % check if angle is larger than 180
    res = beta >= 180;

end

function [startIdx,middleIdx,endIdx] = ...
    aux_incrementIndices(startIdx,middleIdx,endIdx,nrCon)
% auxiliary function for 2D case
% increment indices by 1; additionally, wrap around at the end of the list
     
    startIdx = middleIdx;
    middleIdx = endIdx;
    endIdx = endIdx+1;

    % wrap around if list length is exceeded
    if endIdx > nrCon
        endIdx = 1;
    end

end

function P = aux_oneHalfspace(P,tol)
% polytope with one halfspace (inequality constraints)

% set properties
P.minHRep.val = true;
P.fullDim.val = true;

if all(P.A == 0) && P.b ~= 0
    % constraint infeasible -> empty set
    P.emptySet.val = true;
    P.bounded.val = true;
else
    % not empty
    P.emptySet.val = false;
    P.bounded.val = false;
end

end

function P = aux_twoHalfspaces(P)
% polytope with two halfspaces (inequality constraints)
% note: halfspaces are normalized to unit length

    % dimension
    n = dim(P);
    
    % compute dot product
    dotprod = P.A(1,:) * P.A(2,:)';

    % parallel or anti-parallel?
    if withinTol(dotprod,1)
        % remove the one with the shorter offset
        [P.b,minIdx] = min(P.b);
        P.A = P.A(minIdx,:);

    elseif withinTol(dotprod,-1) && -P.b(2) > P.b(1)
        % -> not feasible
        P = polytope.empty(n);

    else
        % -> feasible
        P.emptySet.val = false;

    end

end

function P = aux_1D(P,type)
% special method for 1D polytopes

    % normalize with respect to b
    P = normalizeConstraints(P,'b');

    % read out properties
    A = P.A; b = P.b;
    Ae = P.Ae; be = P.be;

    if ~isempty(Ae) && (strcmp(type,'all') || strcmp(type,'Ae'))
        % normalization yields equality constraints of the form
        %   ax = 0, ax = 1

        % more than one unique equality constraint is not feasible
        if size(Ae,1) > 1
            if any(withinTol(be,1)) && any(withinTol(be,0))
                % any combination Ax = 1 and Ax = 0 -> empty
                P = polytope.empty(1); return
            end
            
            % indices where b = 1 and b = 0
            idx_1 = withinTol(be,1);
            Ae_1 = Ae(idx_1);
            idx_0 = withinTol(be,0);
            Ae_0 = Ae(idx_0);
    
            if nnz(idx_1) > 1 && ~all(withinTol(Ae_1,Ae_1(1)))
                % more than one constraint with ... = 1
                % and not the same value in A -> empty
                P = polytope.empty(1);
                return
            elseif nnz(idx_0) > 1 && ~all(withinTol(Ae_0,Ae_0(1)))
                % more than one constraint with ... = 0
                % and not the same value in A -> empty
                P = polytope.empty(1);
                return
            end

            % all equality constraints are the same as the first one; take
            % this one and normalize w.r.t Ae
            Ae = 1; be = be(1) / Ae(1);
        end
    end

    if ~isempty(P.A) && (strcmp(type,'all') || strcmp(type,'A'))
        % normalize rows of A matrix
        P = normalizeConstraints(P,'A');
        A = P.A; b = P.b;
        be_ = P.be;
    
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
            P = polytope.empty(1);
            return
        elseif ~isempty(be)
            % there are equality constraints
            if any(be_ > b(1)) || any(be_ < -b(2))
                % additionally, an equality constraint that does not comply
                % with the inequality constraints
                P = polytope.empty(1);
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
%     if ~isempty(A) && ~isempty(Ae)
%         % all are normalized w.r.t to A, so check b values
%         if A(1) == 1 && b(1) 
% 
%     end

    % construct resulting polytope
    P.A = A; P.b = b; P.Ae = Ae; P.be = be;
end

function P = aux_removeAlignedIneq(P)
% only goes through inequality constraints

    % get object properties
    A = P.A;
    b = P.b;
    
    % sort the marix rows to detect aligned normal vectors
    [A,ind] = sortrows(A);
    b = b(ind);
    
    % remove all aligned halfspaces
    counter = 1;
    cTemp = 1;
    A_ = zeros(size(A));
    b_ = zeros(size(b));
    
    while counter < size(A,1)
       
        % check if two normal vectors are identical
        if sum(abs(A(counter,:)-A(counter+1,:))) < eps
            
            a_ = A(counter,:);
            
            % determine minimum b
            bmin = min(b(counter),b(counter+1));
            counter = counter + 2;
            
            while counter <= size(A,1)
                
                if sum(abs(A(counter,:)-a_)) < eps
                   bmin = min(bmin,b(counter));
                   counter = counter + 1;
                else
                   break; 
                end
            end
            
            % store unified normal vector
            A_(cTemp,:) = a_;
            b_(cTemp,:) = bmin;
            cTemp = cTemp + 1;
            
        else
            A_(cTemp,:) = A(counter,:);
            b_(cTemp,:) = b(counter,:);
            cTemp = cTemp + 1;
            counter = counter + 1;
        end  
    end
    
    % add last element
    if counter == size(A,1)
        A_(cTemp,:) = A(end,:);
        b_(cTemp,:) = b(end);
        cTemp = cTemp + 1;
    end
    
    % override constraint set
    P.A = A_(1:cTemp-1,:);
    P.b = b_(1:cTemp-1);
end

% ------------------------------ END OF CODE ------------------------------
