function V = lcon2vert(A,b,varargin)
% lcon2vert - compute the vertices of a polytope
%
% Syntax:
%    V = lcon2vert(A,b,Aeq,beq,TOL,checkbounds)
%
% Inputs:
%    A - inequality constraint matrix ( A*x <= b )
%    b - inequality constraint vector ( A*x <= b )
%    Aeq - (optional) equality constraint matrix ( Aeq*x == beq )
%    beq - (optional) equality constraint vector ( Aeq*x == beq )
%    TOL - (optional) tolerance for vertex computation, default: 1e-10
%    checkbounds - (optional) check if the polytope is bounded before the
%                  vertices are computed (true/false, default = true)
%
% Outputs:
%    V - matrix containing the vertices as column vectors
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Michael Kleder, Matt Jacobson, Niklas Kochdumper
% Written:       09-May-2018
% Last update:   08-December-2022 (MW, polish syntax)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % set default values
    [Aeq,beq,TOL,checkbounds] = setDefaultValues({[],[],1e-10,true},varargin);

    % check if this can be removed:
    if isempty(TOL)
        TOL = 1e-10;
    end
  
    % check number of input arguments
    if nargin == 0
        throw(CORAerror('CORA:notEnoughInputArgs',1));
    elseif nargin == 3
        throw(CORAerror('CORA:specialError',...
            'Since argument Aeq specified, beq must also be specified'));
    elseif nargin > 6
        throw(CORAerror('CORA:tooManyInputArgs',6));
    end
  
    % flatten vectors
    b = b(:); beq = beq(:);
    
    % inputs regarding constraints have to come in pairs
    if xor(isempty(A), isempty(b))
        throw(CORAerror('CORA:specialError',...
            'Since argument A specified, b must also be specified'));
    elseif xor(isempty(Aeq), isempty(beq))
        throw(CORAerror('CORA:specialError',...
            'Since argument Aeq specified, beq must also be specified'));
    end
    
    % ensure that dimension of inequality and equality constraints is equal
    if ~isempty(A) && ~isempty(Aeq) && (size(A,2) ~= size(Aeq,2))
        % dimension of inequality and equality constraints do not match
        throw(CORAerror('CORA:dimensionMismatch',A,Aeq));
    end
    
    % normalize constraints such that norm of each row in A is 0 or 1
    [A,b] = aux_rownormalize(A,b);
    [Aeq,beq] = aux_rownormalize(Aeq,beq);
    
    % extract all inequality constraints that describe equality constraints
    [A,b,Aeq,beq] = aux_extractEqConstr(A,b,Aeq,beq,TOL);
    
    % check if inequality/equality constraints are given
    inequalityConstrained = ~isempty(A);  
    equalityConstrained = ~isempty(Aeq);
    
    % solve equality constraints        
    if equalityConstrained
    
        % compute null space
        Neq = null(Aeq);
        % compute a solution to the set of equality constraints
        x0 = pinv(Aeq)*beq;
    
        if norm(Aeq*x0 - beq) > TOL*norm(beq)
            % infeasible
            V = []; return
    
        elseif isempty(Neq)
    
            Ax0 = A*x0;
            if inequalityConstrained && ~all( Ax0 < b | withinTol(Ax0,b,TOL) )
                % solution to set of equality constraints does not satisfy
                % the inequality constraints -> infeasible
                V = []; return
            else
                % inequality constraints are all satisfied, including
                % vacuously... transpose
                V = x0(:).'; return
            end
        end
    end
    
    
    if inequalityConstrained && equalityConstrained
        % transform inequality constraints to the null-space from equality
        % constraints
        b = b - A*x0;
        A = A*Neq;
    
    elseif equalityConstrained && ~inequalityConstrained
        throw(CORAerror('CORA:specialError',...
            ['Non-bounding constraints detected. '...
            '(Consider box constraints on variables.)']));

    end
    
    % remove redundant constraints
    P = polytope(A,b);
    P = compact_(P,'all',1e-9);
    % read out inequality constraints
    A = P.A;
    b = P.b;
  
    % calculate the vertices from the inequality constraints
    if size(A,2) == 1
        % Special case: 1D
        
        idx_u = sign(A) == 1;
        idx_l = sign(A) == -1;
        idx_0 = sign(A) == 0;
        
        Q = b ./ A;
        U = Q; 
        U(~idx_u) = Inf;
        L = Q;
        L(~idx_l) = -Inf;
        
        [ub,uloc] = min(U);
        [lb,lloc] = max(L);
     
        if ~all(b(idx_0) > 0 | withinTol(b(idx_0),0)) || ub < lb
            % infeasible
            V=[]; return
            % previously: nr=[]; nre=[];
         
        elseif ~isfinite(ub) || ~isfinite(lb)      
            throw(CORAerror('CORA:specialError',...
                'Non-bounding constraints detected. (Consider box constraints on variables.)'));
        end
      
        Zt = [lb;ub];
        
        % previously:
%         if nargout > 1
%             nr=unique([lloc,uloc]); nr=nr(:);
%         end
    else
        % Convert inequality constraints   
        Zt = aux_con2vert(A,b,TOL,checkbounds); 
    end
    
    % convert vertices from null-space back to original space
    if equalityConstrained && ~isempty(Zt)     
        V = bsxfun(@plus, Zt*Neq.', x0(:).');      
    else     
        V = Zt;  
    end
  
end
  
  
% Auxiliary functions -----------------------------------------------------

function [V,nr] = aux_con2vert(A,b,TOL,checkbounds)
% con2vert - convert a convex set of constraint inequalities into the set
%     of vertices at the intersections of those inequalities; i.e., solve
%     the "vertex enumeration" problem. Additionally, identify redundant
%     entries in the list of inequalities.
% 
% Syntax:
%     V = aux_con2vert(A,b)
%     [V,nr] = aux_con2vert(A,b)
% 
% Inputs:
%     A - constraint matrix (m x n with m constraints and n dimensions)
%     b - constraint vector (m x 1 with m constraints)
%
% Outputs:
%     V - vertices (p x n, where p is the number of vertices and n is the
%         dimension)
%     nr - list of the rows in A which are NOT redundant constraints
% 
% NOTES: (1) This program employs a primal-dual polytope method.
%        (2) In dimensions higher than 2, redundant vertices can
%            appear using this method. This program detects redundancies
%            at up to 6 digits of precision, then returns the
%            unique vertices.
%        (3) Non-bounding constraints give erroneous results; therefore,
%            the program detects non-bounding constraints and returns
%            an error. You may wish to implement large "box" constraints
%            on your variables if you need to induce bounding. For example,
%            if x is a person's height in feet, the box constraint
%            -1 <= x <= 1000 would be a reasonable choice to induce
%            boundedness, since no possible solution for x would be
%            prohibited by the bounding box.
%        (4) This program requires that the feasible region have some
%            finite extent in all dimensions. For example, the feasible
%            region cannot be a line segment in 2-D space, or a plane
%            in 3-D space.
%        (5) At least two dimensions are required.
%        (6) See companion function VERT2CON.
%        (7) ver 1.0: initial version, June 2005
%        (8) ver 1.1: enhanced redundancy checks, July 2005
%        (9) Written by Michael Kleder,
%            Modified by Matt Jacobson - March 30, 2011


    % Check if the polytope is bounded
    if checkbounds
        
        [~,bb,~,bbeq] = aux_vert2lcon(A,TOL);
    
        if any(bb<=0) || ~isempty(bbeq)
            throw(CORAerror('CORA:specialError',...
                'Non-bounding constraints detected. (Consider box constraints on variables.)'));
        end
    
        clear bb bbeq
    end

   
    % Initialization ------------------------------------------------------
    
    % determine a point that is located inside the polytope and that is far
    % enough away from the facets of the polytope
    
    % dimension in which the polytope is defined
    n = size(A,2);
    
    if aux_pointInPoly_interior(b,TOL)  
       
        c = zeros(n,1);
   
    else
            
        % slackfun > 0 => point c is located inside the polytope
        % slackfun = 0 (at any entry) => point c is located on the boundary
        slackfun = @(x) b - A*x;

        % Initializer 0
        c = pinv(A)*b;
        s = slackfun(c);

        % use three different methods to compute a point inside the
        % polytope; take 'best' one afterward to avoid refine-branch below
        if ~aux_pointInPoly_approx(s,TOL)
            c = [aux_Initializer1(A,b,c,TOL),...
                 aux_Initializer2(A,b,c,TOL),...
                 center(polytope(A,b))];
            s = slackfun(c);

            % no point inside the polytope could be found
            if ~any(aux_pointInPoly_approx(s,TOL))
                throw(CORAerror('CORA:specialError',...
                    ['Unable to locate a point near '...
                    'the interior of the feasible region.']));
            end

            % choose best one via largest minimum value in s
            [~,idx] = max(min(s,[],1));
            s = s(:,idx);
            c = c(:,idx);
            
        end
        
        % Refinement: If the determined point is located too close to the
        % polytope surface, push it to the interior to increase the
        % numerical stability
        refine = ~aux_pointInPoly_interior(s,TOL);
        if refine

            % determine the surfaces that are too close to the point c
            TOL_ = max(s)*TOL;
            idx = abs(s) < TOL_ | withinTol(abs(s),0,TOL_);

            % Treat the close surfaces as equality constraints
            Amod = A; bmod = b; 
            Amod(idx,:) = []; 
            bmod(idx) = [];

            Aeq = A(idx,:);
            beq = b(idx);

            % Calculate the vertices to the newly generated polytope which
            % is a subspace of the original polytope
            if all(size(A) == size(Amod)) && all(withinTol(A,Amod))
                throw(CORAerror('CORA:specialError',...
                    'Could not find face vertices.'));
            end
            faceVertices = lcon2vert(Amod,bmod,Aeq,beq,TOL,1);
            
            if isempty(faceVertices)
                throw(CORAerror('CORA:specialError',...
                    'Could not find face vertices. Possibly polyhedron is unbounded.'));
            end

            % loop over all possible vertices (choose a feasible vertex)
            foundSolution = false;

            for i = 1:size(faceVertices,1)
                
                try
                    
                    % find the local recession cone vector 
                    c = faceVertices(i,:).';
                    s = slackfun(c);
                    
                    TOL_ = max(s)*TOL;
                    idx = abs(s) < 0 | withinTol(abs(s),0,TOL_);
                    
                    Asub = A(idx,:);
                    bsub = b(idx,:);
                    
                    [aa,bb,aaeq,bbeq] = aux_vert2lcon(Asub);
                    aa = [aa;aaeq;-aaeq];
                    bb = [bb;bbeq;-bbeq];
                    
                    clear aaeq bbeq
                    
                    [bmin,idx] = min(bb);
                    
                    if bmin > -TOL || withinTol(bmin,0,TOL)
                        throw(CORAerror('CORA:specialError',...
                            'We should have found a recession vector (bb<0).'));
                    end
                    
                    
                    % find intersection of polytope with line through facet
                    % centroid
                    Aeq2 = null(aa(idx,:)).';
                    beq2 = Aeq2*c;  
                    

                    linetips = lcon2vert(A,b,Aeq2,beq2,TOL,true);
                    
                    if size(linetips,1) < 2
                        throw(CORAerror('CORA:specialError',...
                            ['Failed to identify line segment through interior. '...
                            'Possibly {x: Aeq*x=beq} has weak intersection with interior({x: Ax<=b}).']));
                    end
                    
                    % take midpoint to the line as the refined point
                    lineCentroid = mean(linetips);
                    
                    clear aa bb
                    
                    c = lineCentroid(:);
                    s = slackfun(c);
                    
                    % suitable point was found -> end search
                    foundSolution = true;
                    break;
                catch ME
                    % no feasible solution found (?)
                end
            end
            
            if ~foundSolution
                throw(CORAerror('CORA:specialError',...
                    'Could not determine a point inside the polytope!'));
            end
        end

        b = s;
    end

   
    % Calculate Vertices --------------------------------------------------
    
    % Normalize the constraints (unit vector b)
    D = bsxfun(@rdivide,A,b); 
    
    % Determine the combinations of inequalies that form a vertex
    try
        k = convhulln(D);
    catch ME
        k = convhulln(D,{'Qs'});
    end
    nr = unique(k(:));
    
    % For each vertex, combine the inequalities that form the vertex and
    % solve a system of linear equations to determine the point
    G = zeros(size(k,1),n);
    ee = ones(size(k,2),1);
    discard = false(1,size(k,1));
    
    for ix = 1:size(k,1)
    
        F = D(k(ix,:),:);
        % check if the inequalities are linearly dependent -> remove
        if aux_lindep(F,TOL) < n
            discard(ix) = true; continue
        end
        
        G(ix,:) = F\ee;
    end
    
    G(discard,:) = [];
    
    % shift the vertices by the initialization point (to original space)
    V = bsxfun(@plus, G, c.'); 
    
    % discard all vertices that are identical
    [~,I] = unique(round(V*1e6),'rows');
    V = V(I,:);
    
end

function [c,fval] = aux_Initializer1(A,b,c,TOL,maxIter)
% Try to determine a point inside the polytope by iterative minimization of
% objective function

    threshold = -10*max(eps(b));
    if nargin > 4
        [c,fval] = fminsearch(@(x) max([threshold;A*x-b]),c,...
            optimset('MaxIter',maxIter));
    else
        [c,fval] = fminsearch(@(x) max([threshold;A*x-b]),c); 
    end
    
end


function c = aux_Initializer2(A,b,c,TOL)
% Determine a point inside the polytope by iterative minimization of
% objective function
    
    % maximum number of iterations
    maxIter = 10000;
    % number of constraints
    nrCon = size(A,1);
    
    % compute psuedo-inverse and augmented constraint matrix
    Ap = pinv(A);        
    Aaug = speye(nrCon) - A*Ap;
    % transpose augmented constraint matrix
    Aaugt = Aaug.';
    
    M = Aaugt*Aaug;
    C = sum(abs(M),2);
    C(C<=0) = min(C(C>0));
    
    slack = b - A*c;
    slack(slack<0) = 0;
     
    IterThresh = maxIter; 
    s = slack; 
    ii = 0;
    
    while ii <= 2*maxIter 
        
        ii = ii + 1; 
        if ii > IterThresh
            IterThresh = IterThresh + maxIter;
        end          
        
        s = s - Aaugt*(Aaug*(s-b))./C;   
        s(s<0) = 0;
        c = Ap*(b-s);
    end
   
end

function [r,idx,Xsub] = aux_lindep(X,varargin)
% Extract a linearly independent set of columns of a given matrix X
%
% Syntax:
%    [r,idx,Xsub] = aux_lindep(X)
%    [r,idx,Xsub] = aux_lindep(X,tol)
%
% Inputs:
%    X - matrix
%    tol - rank estimation tolerance (default = 1e-10)
%
% Outputs:
%    r - rank estimate
%    idx - indices (into X) of linearly independent columns
%    Xsub - extracted linearly independent columns of X

    % set default values
    tol = setDefaultValues({1e-10},varargin{:});

    if ~nnz(X)
        % X has no non-zeros and hence no independent columns
        Xsub = []; idx = []; return
    end    
    
    % QR decomposition
    [~,R,E] = qr(X,0);
    diagr = abs(diag(R));
    
    % Rank estimation
    r = find(diagr >= tol*diagr(1),1,'last');
    
    % set optional output arguments
    if nargout > 1
        idx=sort(E(1:r));
        idx=idx(:);
        if nargout > 2
            Xsub=X(:,idx);
        end
    end                

end
     
function [A,b] = aux_rownormalize(A,b)
% Modifies A,b data pair so that norm of each row in A is either 0 or 1

    % return if no constraints given
    if isempty(A)
        return
    end
    
    normsA = sqrt(sum(A.^2,2));
    idx = normsA>0;
    A(idx,:) = bsxfun(@rdivide,A(idx,:),normsA(idx));
    b(idx) = b(idx) ./ normsA(idx);

end

function res = aux_pointInPoly_approx(s,TOL)
% Determines if a point is approximately (up to tolerance) located inside 
%    the polytope; since
%        s(x) = b - A*x,  with s \in R^m (m ... number of constraints) 
%    and we need to fulfill
%        A*x <= b + tol <=>  -tol <= b - A*x  <=>  -tol <= s(x)
%    the point is contained if all constraints are fulfilled, i.e.,
%        \forall i \in {1,...,m}: -tol <= s_i(x)
%    or definitely not contained if
%        \exists i \in {1,...,m}: -tol > s_i(x)
% 
% Syntax:
%    res = approxinpoly(s,TOL)
%
% Inputs:
%    s - slack value of the point x (= b - A*x from inequality A*x <= b)
%    TOL - tolerance
%
% Outputs:
%    res - true/false
    
    if any(-TOL > s)
        % one violation found
        res = false;
    else
        % adapt tolerance
        TOL_ = abs(max(max(s))*TOL);
        % check all inequalities
        res = all( 0 < s | withinTol(s,0,TOL_) );
    end

    % old version:
%     smax = max(s);   
%     if smax <= 0
%         res = false;
%     else
%         res = all(s >= -smax*TOL);
%     end
    
end
   
function res = aux_pointInPoly_interior(s,TOL)
% Determines if a point is an interior point of the polytope (up to
%    tolerance); since
%        s(x) = b - A*x,  with s \in R^m (m ... number of constraints) 
%    and we need to fulfill
%        A*x < b + tol <=>  -tol < b - A*x  <=>  -tol < s(x)
%    the point is contained if all constraints are fulfilled, i.e.,
%        \forall i \in {1,...,m}: -tol < s_i(x)
%    or definitely not contained if
%        \exists i \in {1,...,m}: -tol >= s_i(x)
% 
% Syntax:
%    res = strictinpoly(s,TOL)
%
% Inputs:
%    s - slack value of the point x (b - A*x from inequality A*x <= b)
%    TOL - tolerance
%
% Outputs:
%    res - true/false

    if any(-TOL > s | withinTol(0,s,TOL))
        res = false;
    else
        res = all(-TOL < s);
    end

    % old version:
%     smax = max(s);
%     if smax <= 0
%         res = false;
%     else
%         res = all(s >= smax*TOL);
%     end

end
   
function [A,b,Aeq,beq] = aux_vert2lcon(V,varargin)
% vert2lcon - an extension of Michael Kleder's vert2con function, used for
%    finding the linear constraints defining a polytope in R^n given its
%    vertices. This wrapper extends the capabilities of vert2con to also
%    handle cases where the polytope is not solid in R^n, i.e., where the
%    polytope is defined by both equality and inequality constraints.
%    Any point x inside the polytope will/must satisfy
%       A*x  <= b
%       Aeq*x = beq
%    up to machine precision issues.
%
% Syntax:
%    [A,b,Aeq,beq] = aux_vert2lcon(V,TOL)
%
% Inputs:
%    V - vertices (N x n, where N is the number of vertices and n the
%                  dimension)
%    tol - rank estimation tolerance (default = 1e-10)
%
% Outputs:
%    A - inequality constraint matrix
%    b - inequality constraint vector
%    Aeq - equality constraint matrix
%    beq - equality constraint vector
%
% Example:
%    % 3D region defined by x+y+z = 1, x >= 0, y >= 0, z >= 0
%    V = eye(3);
%    [A,b,Aeq,beq] = aux_vert2lcon(V);
%
%    % output
%    A =
%         0.4082   -0.8165    0.4082
%         0.4082    0.4082   -0.8165
%        -0.8165    0.4082    0.4082
%    b =
%         0.4082
%         0.4082
%         0.4082
%    Aeq =
%         0.5774    0.5774    0.5774
%    beq =
%         0.5774

    % set default values
    tol = setDefaultValues({1e-10},varargin);
    
    % number of vertices and dimension
    [M,N] = size(V);
    
    if M == 1
        A = []; b = [];
        Aeq = eye(N); beq=V(:);
        return
    end
    
    p = V(1,:).';
    X = bsxfun(@minus,V.',p);
    
    % from now on, we need Q to be full column rank and prefer E to be compact
    
    if M > N
        % X is wide

        % QR decomposition
        [Q, R, E] = qr(X,0);  
        % ...economy-QR ensures that E is compact
        % ...Q automatically full column rank since X wide
    
    else
        % X is tall, hence non-solid polytope
    
        % QR decomposition
        [Q, R, P] = qr(X);
        % ...non-economy-QR so that Q is full-column rank
    
        % No way to get E compact, this is the alternative:
        [~,E] = max(P);
    
        clear P
    end 
    
    diagr = abs(diag(R));
    
    
    if nnz(diagr)    
    
        % rank estimation
        r = find(diagr >= tol*diagr(1), 1, 'last');
        
        iE = 1:length(E);
        iE(E) = iE;
        
        Rsub = R(1:r,iE).';
    
        if r > 1
        
            [A,b] = aux_vert2con(Rsub,tol);
        
        elseif r == 1
        
            A = [1;-1];
            b = [max(Rsub);-min(Rsub)];
    
        end
        
        % inequality constraints
        A = A*Q(:,1:r).';
        b = bsxfun(@plus,b,A*p);
        
        % equality constraints
        if r < N
            Aeq = Q(:,r+1:end).';      
            beq = Aeq*p;
        else
            Aeq = [];
            beq = [];
        end
    
    else
        % rank = 0 -> all points are identical
    
        % no inequality constraints
        A = []; b = [];
        % only equality constraints
        Aeq = eye(N); beq = p;
    end
   
end
           
function [A,b] = aux_vert2con(V,tol)
% VERT2CON - convert a set of points to the set of inequality constraints
%            which most tightly contain the points; i.e., create
%            constraints to bound the convex hull of the given points
%
% Syntax:
%     [A,b] = aux_vert2con(V)
%     [A,b] = aux_vert2con(V,tol)
%
% Inputs:
%     V - a set of points, each ROW of which is one point
%
% Outputs
%     A,b - a set of constraints such that A*x <= b defines the region of
%           of space enclosing the convex hull of the given points
%
% For dimension n:
%     V = p x n matrix (p vertices in dimension n)
%     A = m x n matrix (m constraints in dimension n)
%     b = m x 1 vector (m constraints)
%
% NOTES: (1) In higher dimensions, duplicate constraints can
%            appear. This program detects duplicates at up to 6
%            digits of precision, then returns the unique constraints.
%        (2) See companion function CON2VERT.
%        (3) ver 1.0: initial version, June 2005.
%        (4) ver 1.1: enhanced redundancy checks, July 2005
%        (5) Written by Michael Kleder,
%            Modified by Matt Jacobson - March 29, 2011

    try
        k = convhulln(V);
    catch
        k = convhulln(V,{'Qs'});
    end
    c = mean(V(unique(k),:));
    
    
    V = bsxfun(@minus,V,c);
    A = NaN(size(k,1),size(V,2));
    
    n = size(V,2);
    ee = ones(size(k,2),1);
    rc = 0;
    
    for ix = 1:size(k,1)
        F = V(k(ix,:),:);
        if aux_lindep(F,tol) == n
            rc = rc + 1;
            A(rc,:) = F\ee;
        end
    end
    
    A = A(1:rc,:);
    b = ones(size(A,1),1);
    b = b + A*c';
    
    % eliminate duplicate constraints:
    [A,b] = aux_rownormalize(A,b);
    [~,I] = unique(round([A,b]*1e6),'rows');
    
    A = A(I,:); % NOTE: rounding is NOT done for actual returned results
    b = b(I);

end
         
         
function [A,b,Aeq,beq] = aux_extractEqConstr(A,b,Aeq,beq,tol)
% Extract all inequality constraints of the form
%    a*x <= b  &&  a*x >= b <=> -a*x <= -b
% as they can be described by an equivalent equality constraint of the form
%    a*x = b
    
    % Heuristic: sort the constraints according to a hash function to make
    % it easier to identify similar constraints
    A_ = abs(A);
    temp = 1:size(A,2);
    [A_,ind] = sortrows([A_*temp',A_]);
    A = A(ind,:);
    b = b(ind,:);
    
    
    % loop over all constraints
    i = 1;
    while i < size(A,1)
        
        % get index of last row that has the same hash value as the current 
        % row
        lInd = i;
        for k = i+1:size(A,1)
            if A_(k,1) ~= A_(i,1)
                break; 
            else
                lInd = lInd + 1;
            end
        end
        
        % compare all rows that have the same hash value to detect possible
        % equality constraints
        indRem = [];
        
        for j = i:lInd
            for k = j+1:lInd
%                 if all(withinTol(abs(A(j,:)+A(k,:)),0,tol)) ...
%                         && withinTol(abs(b(j)+b(k)),0,tol)
                if all(abs(A(j,:)+A(k,:)) < tol) && abs(b(j)+b(k)) < tol
                    indRem = [indRem;j;k];
                    Aeq = [Aeq;A(j,:)];
                    beq = [beq;b(j)];
                end
                cond_old = all(abs(A(j,:)+A(k,:)) < tol) && abs(b(j)+b(k)) < tol;
                cond_new = all(withinTol(abs(A(j,:)+A(k,:)),0,tol)) && withinTol(abs(b(j)+b(k)),0,tol);
                if cond_old ~= cond_new
                    disp(" check methods ");
                end
            end
        end
        
        indRem = unique(indRem);
        A(indRem,:) = [];
        b(indRem) = [];
        A_(indRem,:) = [];
        
        i = lInd + 1 - length(indRem);
    end
end

% ------------------------------ END OF CODE ------------------------------
