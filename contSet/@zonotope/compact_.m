function Z = compact_(Z,method,tol,varargin)
% compact_ - returns equal zonotope in minimal representation
%
% Syntax:
%    Z = compact_(Z,method,tol)
%
% Inputs:
%    Z - zonotope object
%    method - method for redundancy removal
%             'zeros': removes all generators from a zonotope
%                  with zeros in all dimensions, so that every generator
%                  the resulting zonotope has at least one non-zero entry
%             'all': combines aligned generators to a single generator;
%                  a tolerance is used to determine alignment, so this
%                  function does not necessarily return an
%                  over-approximation of the original zonotope---for this,
%                  use zonotope/reduce instead
%    tol - tolerance
%
% Outputs:
%    Z - zonotope object
%
% Example:
%    Z1 = zonotope([0;0],[1 0 -2 0 3 4; 0 0 1 0 -2 1]);
%    Z2 = compact(Z1);
%    
%    plot(Z1); hold on;
%    plot(Z2,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/compact, zonotope/reduce

% Authors:       Mark Wetzlinger, Matthias Althoff
% Written:       15-January-2009
% Last update:   27-August-2019
%                05-October-2024 (MW, remove superfluous 'aligned', rewrite deleteAligned)
%                27-May-2025 (TL, bug fix, outputs are now consistent for 'all')
% Last revision: 29-July-2023 (MW, merged from deleteZeros/deleteAligned)

% ------------------------------ BEGIN CODE -------------------------------

switch method
    case 'zeros'
        Z = aux_deleteZeros(Z,tol);
    case 'all'
        Z = aux_deleteAligned(Z,tol);
end

end


% Auxiliary functions -----------------------------------------------------

function Z = aux_deleteZeros(Z,tol)
    
    % filter zero generators
    Z.G = nonzeroFilter(Z.G,tol);

end

function Z = aux_deleteAligned(Z,tol)
    
    % delete zero-generators
    G = nonzeroFilter(Z.G,tol);

    % quick exit for 1D case
    if dim(Z) == 1
        Z.G = sum(abs(G));
        return
    end
    
    % normalize generators
    vnormG = vecnorm(G);
    G_norm = G ./ vnormG;
    nrGen = size(G,2);

    % sort generators for consistent outputs
    % TL: We had issues with it compacting less generators depending on the
    % order of the generators. Sorting them should solve the issue and
    % should not be slower than the subsequent computations.
    [~,idx] = sort(vnormG);
    G_norm = G_norm(:,idx);
    G = G(:,idx);
    
    % find aligned generators: since the generators are normalized, aligned
    % generators must have a dot product of 1 (parallel) or -1
    % (anti-parallel) with the generator in question; generators are then
    % unified by addition (taking 1/-1 into account)
    idxKeep = true(1,nrGen);
    for i=1:nrGen
        % only check those that have not yet been removed in a previous
        % iteration
        if idxKeep(i)
            % compute dot product
            dotprod = G_norm(:,i)' * G_norm;

            % check for which pairs the value is 1 or -1
            idxParallel = withinTol(dotprod,1,tol);
            idxAntiParallel = withinTol(dotprod,-1,tol);

            % add all generators with dotprod = 1 to ith generator,
            % subtract all generators with dotprod = -1 from ith generator
            % (note: dotprod with itself is 1, therefore no extra "+G(:,i)")
            if nnz(idxParallel) > 1 || any(idxAntiParallel)
                G(:,i) = sum(G(:,idxParallel),2) ...
                        - sum(G(:,idxAntiParallel),2);

                % ith generator and all (anti-)parallel generators have been
                % subsumed -> remove them from further checks; keep ith
                % generator for final removal after loop
                idxParallel(i) = false;
                idxKeep = idxKeep & ~idxParallel & ~idxAntiParallel;
            end
        end
    end

    % lengthening has already taken place in the loop above, so we only
    % need to remove all subsumed generators
    Z.G = G(:,idxKeep);

end

% ------------------------------ END OF CODE ------------------------------
