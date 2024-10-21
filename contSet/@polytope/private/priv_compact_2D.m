function [A,b,Ae,be,empty] = priv_compact_2D(A,b,Ae,be,tol)
% priv_compact_2D - special computation for 2D case without linear
%    programs: First, we order the (already normalized) constraints by
%    angle. Then, we check whether the middle constraint between two
%    adjacent constraints is required by computing the intersecting vertex
%    between the outer two constraints and checking whether that middle
%    constraint has a smaller support function value than the intersecting
%    vertex in the direction of the middle constraint
%
% Syntax:
%    [A,b,Ae,be,empty] = priv_compact_2D(A,b,Ae,be,tol)
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
    % remove duplicates
    [Ae,be,empty] = priv_compact_alignedEq(Ae,be,tol);
    % emptiness already determined
    if empty
        A = []; b = []; Ae = []; be = [];
        return
    end
end


if ~isempty(A)

    % transpose constraints, read out number
    A = A';
    nrCon = length(b);

    % constraints are normalized, order by angle
    angles = atan2d(A(2,:),A(1,:));
    [angles,idx] = sort(angles);
    A = A(:,idx);
    b = b(idx);

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
    if length(b) == 1
        A = A';
        return
    elseif length(b) == 2
        [A,b,empty] = aux_twoHalfspaces(A',b);
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
    A = A(:,idxKept)';
    b = b(idxKept);

    % TODO: check if empty without linprog (?)
end

end


% Auxiliary functions -----------------------------------------------------

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

function [A,b,empty] = aux_twoHalfspaces(A,b)
% polytope with two halfspaces (inequality constraints)
% note: halfspaces are normalized to unit length

empty = false;

% compute dot product
dotprod = A(1,:) * A(2,:)';

% parallel or anti-parallel?
if withinTol(dotprod,1)
    % remove the one with the shorter offset
    [b,minIdx] = min(b);
    A = A(minIdx,:);

elseif withinTol(dotprod,-1) && -b(2) > b(1)
    empty = true;

end

end

% ------------------------------ END OF CODE ------------------------------
