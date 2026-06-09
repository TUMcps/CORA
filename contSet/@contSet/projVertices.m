function V = projVertices(S,varargin)
% projVertices - computes the vertices of a 2D projection of a set based on
%    support function evaluation of the projected set; if no more support
%    vectors can be found, the algorithm terminates
%    this function also supports degenerate sets (lines, points)
%
% Syntax:
%    V = projVertices(S)
%    V = projVertices(S,dims)
%    V = projVertices(S,dims,alg)
%
% Inputs:
%    S - contSet object
%    dims - dimensions for projection
%    alg - algorithm, which is either "angle" (default) or "supportFunc" 
%
% Outputs:
%    V - list of vertices in the projected space
%
% Example:
%    Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
%    A = [1 1 1]; b = 1;
%    cZ = conZonotope(Z,A,b);
%
%    V = projVertices(cZ);
%
%    figure; hold on;
%    plot(cZ);
%    scatter(V(1,:),V(2,:),16,'r','filled');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger, Niklas Kochdumper
% Written:       21-December-2022
% Last update:   29-April-2024 (TL, increased tol for init duplicate check)
%                20-October-2025 (NK, added new faster algorithm)
%                19-February-2025 (TL, idx shift bug fix)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % too many input arguments
    narginchk(1,3);
    
    % set default values
    [dims,alg] = setDefaultValues({[1,2],'angle'},varargin);
    
    % check input arguments (only supported for convex sets with edges)
    inputArgsCheck({ {S,'att',{'conZonotope','interval','zonoBundle','zonotope','polytope'}} ...
                     {dims,'att','numeric',{'integer','positive',@(dims) numel(dims) == 2}} ...
                     {alg,'str',{'angle','supportFunc'}}});

    % compute the vertices of the projection to 2D (since the faster 
    % "angle" algorithm might fail for special cases like degenerate sets, 
    % etc., we call the more stable "supportFunc" algorithm in case of 
    % failure)
    if strcmp(alg,'angle')
        try
            V = aux_projVerticesAngle(S,dims);
        catch
            V = aux_projVerticesSupportFunc(S,dims);
        end
    else
        V = aux_projVerticesSupportFunc(S,dims);
    end
end


% Auxiliary functions -----------------------------------------------------

function V = aux_projVerticesAngle(S,dims)
% computes the vertices of a 2D projection of a convex set based on the 
% angle of edged. The algorithm works similar like the one for plotting
% zonotopes, where starting from one extrem point we walk around the convex
% set along the edges, where the best edge to walk on is
% determined based on the angle.

    tol = 1e-10;

    % if polytope vertex representation given, much simpler projection
    if isa(S,'polytope') && S.isVRep.val
        V = S.V(dims,:);
        ind = convhull(V(1,:),V(2,:));
        V = V(:,ind(1:end-1));
        return
    end

    % suppress warning for close to singular matrices in equation systems
    w = warning(); warning('off');

    % convert zonotope bundles to constrained zonotopes
    if isa(S,'zonoBundle') || isa(S,'zonotope') || isa(S,'interval')
        S = conZonotope(S);
    end

    % center set at the origin
    c = center(S);
    S = S + (-c);

    % extract constraints for different types of set representations
    % (common form of sets:  S = c_x + G_x*a,  A*a <= b,  Aeq*a <= beq)
    if isa(S,'conZonotope')
        c_x = S.c; G_x = S.G;
        m = size(S.G,2); A = [eye(m);-eye(m)]; b = ones(2*m,1);
        Aeq = S.A; beq = S.b;
    elseif isa(S,'polytope')
        S = compact(S,'AtoAe');
        n = dim(S); c_x = zeros(n,1); G_x = eye(n);
        A = S.A; b = S.b; Aeq = S.Ae; beq = S.be;
    end

    % project the set to the desired dimensions
    c_x = c_x(dims); G_x = G_x(dims,:);

    % compute null-space of the equality constraints
    if ~isempty(Aeq)
        T = null(Aeq);
        p = pinv(Aeq)*beq;
    else
        T = eye(size(G_x,2)); p = zeros(size(G_x,2),1);
    end

    % construct polytope for inequality constraints in the nullspace    
    P = compact(polytope(A*T,b - A*p),'aligned');
    P_A = P.A; P_b = P.b;

    % determine maximum point in x-direction (and make sure it is a vertex)
    f = G_x*T; f = -f(1,:);
    problem.f = f'; problem.Aineq = P_A; problem.bineq = P_b;
    problem.Aeq = []; problem.beq = []; cnt = 1;

    while true
       
        % compute maximum point using linear programming
        [fac,~,exitflag] = CORAlinprog(problem);
    
        if exitflag < 0
            throw(CORAerror('CORA:solverIssue'));
        end

        % check if point is a vertex and modify optimization problem if not
        dist = abs(P_A*fac - P_b);
        ind = find(dist < tol);

        if (length(ind) >= size(T,2) && rank(P_A(ind,:)) == size(T,2)) ...
                                                            || cnt > 10
            break;
        else
            problem.Aeq = [problem.Aeq;problem.Aineq(ind,:)];
            problem.beq = [problem.beq;problem.bineq(ind,:)];
            problem.Aineq(ind,:) = []; problem.bineq(ind,:) = [];
            problem.f = problem.f + (1:length(problem.f))';
        end

        cnt = cnt + 1;
    end

    x = c_x + G_x*(T*fac + p);

    % loop until a full circle is complete
    V = x; opts.UT = true; q = length(fac); cnt = 1; checkEnd = false;

    while true

        % catch case involving a degenerate set
        if norm(x) < tol
            throw(CORAerror("CORA:specialError", ...
                                    "Degenerate sets are not supported."));
        end

        % compute angle of current point together with rotation matrix
        phi_x = atan2(x(2),x(1));
        Rot = [cos(phi_x) -sin(phi_x); sin(phi_x) cos(phi_x)];
    
        % find constraints that are active at this point
        dist = abs(P_A*fac - P_b);
        indAll = find(dist < tol);

        % check if a full circle is complete
        if size(V,2) == 1
            indFirst = indAll;
        elseif all(ismember(indAll,indFirst)) || ...
               all(ismember(indFirst,indAll)) || ...
               (checkEnd && norm(x - V(:,1)) < tol)
            break;
        end

        % loop over all possible combinations of n halfspaces that might 
        % form the current vertex 
        comb = combinator(length(indAll),length(fac),'c');

        if size(comb,1) > 1000
            throw(CORAerror("CORA:specialError", ...
                                          "Computationally infeasible."));
        end

        q_ = q*size(comb,1); valid = zeros(q_,1); cntFac = 1; 
        D = zeros(2,q_); Dfac = zeros(q,q_); phi = zeros(q_,1);
        recomp = false;
        
        for j = 1:size(comb,1)

            ind = indAll(comb(j,:));
    
            % transform equation system into upper-triangular form using
            % LU-decomposition (since it has to be solved multiple times)
            [Q,R] = qr(P_A(ind,:));
    
            % recompute the factor values for the current vertex to avoid
            % amplification of numerical errors
            if ~recomp && all(abs(diag(R)) > 1000*eps)
                fac = linsolve(R,Q'*P_b(ind),opts);
                recomp = true;
            end
    
            % consider all potential combinations of constraints                
            for i = 1:q
        
                % determine direction of the edge            
                b_ = P_b(ind); b_(i) = b_(i) + 1;
                x_ = linsolve(R,Q'*b_,opts);
                d = x_ - fac; d = d/norm(d);
    
                if max(P_A(ind,:)*(fac + d) - P_b(ind)) > ...
                                      max(P_A(ind,:)*(fac - d) - P_b(ind))
                    d = -d;
                end
    
                Dfac(:,cntFac) = d;

                % check if the edge is really an edge of the polytope
                % (required since there can be redundant halfspaces)
                if max(P_A(indAll,:)*(fac + d) - P_b(indAll)) < 1e-10
                    valid(cntFac) = 1; 
                end
    
                % transform edge from factor space to the state space
                d = G_x*T*d;
                D(:,cntFac) = d;

                % exclude 0 vectors in the state space
                if norm(d) < 10*eps
                    valid(cntFac) = 0;
                end
        
                % compute angle in the 2D plane
                d_ = Rot'*d(1:2);
                phi(cntFac) = atan2(d_(2),d_(1));

                cntFac = cntFac + 1;
            end
        end
    
        % determine best edge based on the angle
        ind1 = find(valid);
        ind2 = find(phi(ind1) >= 0);
        [~,ind3] = min(phi(ind1(ind2)));
    
        dfac = Dfac(:,ind1(ind2(ind3)));
    
        % compute next point by walking along the edge until another
        % inequality constraint becomes active
        ind_ = find(dist >= tol);
        lamda = inf*ones(length(ind_),1);
    
        for i = 1:length(ind_)
            a_ = P_A(ind_(i),:); b_ = P_b(ind_(i));
            dirProjection = a_*dfac;
            if abs(dirProjection) > eps
                lamda(i) = (b_ - a_*fac)/dirProjection;
            end
        end
    
        ind1 = find(lamda >= 0);
        [~,ind2] = min(lamda(ind1));
    
        fac = fac + dfac * lamda(ind1(ind2));
        x = c_x + G_x*(T*fac + p);
        V = [V,x];

        % check if algorithm fails to converge
        cnt = cnt + 1;

        if norm(x - V(:,1)) > 100*tol
            checkEnd = true;
        end

        if cnt > 1e5
            throw(CORAerror("CORA:specialError", ...
                                        "Algorithm failed to converge."));
        end
    end

    % shift by center
    V = V + c(dims);
    V = V(1:2,1:end-1);

    warning(w);
end

function V = aux_projVerticesSupportFunc(S,dims)
% computes the vertices of a 2D projection of a set based on support 
% function evaluation of the projected set; if no more support vectors can 
% be found, the algorithm terminates this function also supports degenerate 
% sets (lines, points)

    % convert zonotope boundles before projecting since projection is not
    % exact for zonotope boundles
    if isa(S,'zonoBundle')
        S = conZonotope(S);
    end
    
    % init vertices
    V = [];
    if representsa_(S,'emptySet',1e-10)
        % no vertices if set is empty
        return
    end
    
    % other options for support function evaluation
    otherOptions = {};
    if isa(S,'polyZonotope') || isa(S,'conPolyZono')
        otherOptions = {'interval',8,1e-3};
    end
    
    % compute support vectors of three directions
    
    % 0 degrees
    dir = zeros(dim(S),1); dir(dims) = [1;0];
    [~,V(:,1)] = supportFunc_(S,dir,'upper',otherOptions{:});
    % plot([0,1],[0,0],'k');
    
    % 120 degrees
    angle_pi = 120*pi / 180;
    dir(dims) = [cos(angle_pi) -sin(angle_pi); sin(angle_pi) cos(angle_pi)] * [1;0];
    [~,V(:,2)] = supportFunc_(S,dir,'upper',otherOptions{:});
    
    % 240 degrees
    angle_pi = 240*pi / 180;
    dir(dims) = [cos(angle_pi) -sin(angle_pi); sin(angle_pi) cos(angle_pi)] * [1;0];
    [~,V(:,3)] = supportFunc_(S,dir,'upper',otherOptions{:});
    
    % copy last point for easier indexing (will be deleted later)
    V(:,4) = V(:,1);
    V = V(dims,:);

    % ensure that there are no duplicates (other than 1 = end)
    if compareMatrices(V(:,3),V(:,1:2),1e-12,'subset')
        % vertex (240 degrees) is duplicate
        V(:,3) = [];
    end
    if withinTol(V(:,2),V(:,1))
        V(:,2) = [];
    end
    
    % all sections between already computed vertices have to be investigated
    sections = mat2cell([(1:size(V,2)-1)', (2:size(V,2))'], ...
                                                    ones(size(V,2)-1,1),2);
    
    % split angles between neighboring directions until no new information
    while ~isempty(sections)
        
    
        % analyze first section in the list of sections
        section = sections{1};
        sections = sections(2:end);
    
        % direction = normal vector on halfspace spanned by points on section
        v = V(:,section(2)) - V(:,section(1));
        dir_ = [v(2); -v(1)];
        dir_ = dir_ / vecnorm(dir_);
        
        % compute support vector for the new direction
        dir(dims) = dir_;
        [~,V_new] = supportFunc_(S,dir,'upper',otherOptions{:});
        V_new = V_new(dims);
        
        % compute vectors:
        % - from start vertex to computed vertex
        % - from computed vertex to end vertex
        ptsStartMidEnd = [V_new-V(:,section(1)), V(:,section(2))-V_new];
    
        % check whether sections is completed
        if compareMatrices(V_new,V,1e-6,'subset') ...
                || rank(ptsStartMidEnd,1e-6) < 2
            % new vertex is on a line with start and end points of the
            % current section -> discard vertex, section completed
        else
            % vertices are not on a line
            
            % add vertex to list at index between points from current section
            V = [V(:,1:section(1)) V_new V(:,section(2):end)];
    
            % shift indices of other sections
            for s=1:length(sections)
                sections{s} = sections{s} + 1;
            end
            % add two new sections at the beginning
            sections = [{[section(1) section(1)+1]; ...
                                    [section(2) section(2)+1]}; sections];
    
        end
    
    end
    
    % remove last vertex (was only added for convenience)
    V = V(:,1:end-1);

    % remove all collinear vertices: support function queries return LP
    % vertices in the full-dimensional space, which can project onto edges
    % (not vertices) of the 2D polygon; check all vertices, not just the
    % initial three
    V = aux_removeCollinear(V);
end

function V = aux_removeCollinear(V)
% remove all vertices that are collinear with their circular neighbors

    if size(V,2) <= 2
        return
    end

    toDelete = false(1, size(V,2));
    nV = size(V,2);

    for k = 1:nV
        % circular neighbor indices
        prev = mod(k-2, nV) + 1;
        next = mod(k, nV) + 1;

        % check collinearity with neighbors
        ptsStartMidEnd = [V(:,k)-V(:,next), V(:,prev)-V(:,k)];
        if rank(ptsStartMidEnd, 1e-6) < 2
            toDelete(k) = true;
        end
    end

    V(:, toDelete) = [];
end

% ------------------------------ END OF CODE ------------------------------
