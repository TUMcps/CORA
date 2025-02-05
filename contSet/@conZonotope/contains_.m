function [res,cert,scaling] = contains_(cZ,S,method,tol,maxEval,certToggle,scalingToggle)
% contains_ - determines if a constrained zonotope contains a set or a
%    point
%
% Syntax:
%    res = contains_(cZ,S,method,tol,maxEval,certToggle,scalingToggle)
%
% Inputs:
%    cZ - conZonotope object
%    S - contSet object or single point or matrix of points
%    method - method used for the containment check.
%       The available options are:
%           - 'exact': Checks for containment by either enumerating the
%               vertices of the inbody (if possible), or enumerating the
%               halfspaces of cZ.
%           - 'approx': Checks for containment using the algorithm from
%               Theorem 1 [1] if S is polytopic. Otherwise, a polytopic
%               over-approximation of S is computed.
%               This algorithm has polynomial runtime, and if
%               res==1, containment is assured. Otherwise, S could still be
%               contained in cZ, but this cannot be verified in polynomial
%               time.
%          For the next methods, note that if both certToggle and
%          scalingToggle are set to 'false', then res will be set to
%          'false' automatically, and the algorithms will not be executed.
%          This is because stochastic/optimization-based algorithms can not
%          confirm containment, so res = true can never happen. However, if
%          maxEval is set high enough, and res = false but cert = false,
%          one might conclude that with good probability, containment
%          holds.
%           - 'sampling:primal': Solves the containment stochastically,
%               using a similar algorithm as the Shenmaier vertex sampling
%               from [4].
%           - 'sampling:dual': Solves the containment stochastically, using
%               a similar algorithm as the Shenmaier halfspace sampling
%               from [4].
%    tol - tolerance for the containment check; the higher the
%       tolerance, the more likely it is that points near the boundary of 
%       cZ will be detected as lying in cZ, which can be useful to 
%       counteract errors originating from floating point errors.
%    maxEval - only, if 'sampling:primal' or 'sampling:dual' is used:
%       Number of maximal function evaluations.
%    certToggle - if set to 'true', cert will be computed (see below).
%    scalingToggle - if set to 'true', scaling will be computed (see
%       below).
%
% Outputs:
%    res - true/false
%    cert - returns true iff the result of res could be
%           verified. For example, if res=false and cert=true, S is
%           guaranteed to not be contained in cZ, whereas if res=false and
%           cert=false, nothing can be deduced (S could still be
%           contained in cZ).
%           If res=true, then cert=true.
%           Note that computing this certification may marginally increase
%           the runtime.
%    scaling - returns the smallest number 'scaling', such that
%           scaling*(cZ - center(cZ)) + center(cZ) contains S.
%           For the methods 'approx' and 'approx:st' this is an upper
%           bound, for 'sampling:primal' and 'sampling:dual', this number
%           is a lower bound. Note that computing this scaling factor may 
%           significantly increase the runtime.
%
% Example: 
%    % generate constrained zonotopes
%    Z = [0 2 -2 1;0 1.5 1 -1.5];
%    A = [1 1 1]; b = 1;
%    cZ1 = conZonotope(Z,A,b);
% 
%    Z = [0 2 0 0;0 1 1 0];
%    A = [1 1 -1]; b = 0;
%    cZ2 = conZonotope(Z,A,b);
%
%    Z = [1 2 0 0;1 1 1 0];
%    A = [1 1 -1]; b = 0;
%    cZ3 = conZonotope(Z,A,b);
%
%    % check for containment
%    contains(cZ1,cZ2)
%    contains(cZ1,cZ3)
%
%    % visualization
%    figure; hold on;
%    plot(cZ1,[1,2],'r');
%    plot(cZ2,[1,2],'b');
%
%    figure; hold on;
%    plot(cZ1,[1,2],'r');
%    plot(cZ3,[1,2],'b');
%
% References:
%    [1] Sadraddini et al.: "Linear Encodings for Polytope Containment
%        Problems", CDC 2019
%    [2] JK Scott, DM Raimondo, GR Marseglia, RD Braatz: "Constrained 
%        zonotopes: A new tool for set-based estimation and fault
%        detection", Automatica 69, 126-136
%        detection, Automatica 69, 126-136
%    [3] Kulmburg A., Brkan I., Althoff M.,: "Search-based and Stochastic
%        Solutions to the Zonotope and Ellipsotope Containment Problems",
%        ECC 2024
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/contains, zonotope/contains_

% Authors:       Niklas Kochdumper, Mark Wetzlinger, Matthias Althoff, Adrian Kulmburg, Ivan Brkan
% Written:       14-November-2019 
% Last update:   15-November-2022 (MW, return logical array for points)
%                25-November-2022 (MW, rename 'contains')
%                09-March-2023 (MA, added reference for point enclosure)
%                05-February-2024 (AK, IB, added sampling methods, cert and scaling)
% Last revision: 27-March-2023 (MW, rename contains_)

% ------------------------------ BEGIN CODE -------------------------------

% No matter what, the algorithms will run better if we simplify cZ as much
% as possible
cZ = compact(cZ);

% Check out trivial cases
[cZ_isPoint, p] = representsa(cZ, 'point');
if cZ_isPoint
    if isnumeric(S)
        res = max(max(abs(S-p))) <= tol;
        cert = true;
        if res
            scaling = 0;
        else
            scaling = Inf;
        end

    else % S is not numeric
        [S_isPoint, q] = representsa(S, 'point');
        if S_isPoint
            res = all(q==p);
            cert = true;
            if res
                scaling = 0;
            else
                scaling = Inf;
            end
        end
    end
    return
end
% cZ is not a point

if representsa(cZ, 'emptySet')
    if isnumeric(S)
        if isempty(S)
            res = true;
            scaling = 0;
            cert = true;
        else
            res = false;
            scaling = Inf;
            cert = true;
        end

    elseif representsa(S, 'emptySet')
        res = true;
        scaling = 0;
        cert = true;
        
    else % S is not numeric and not empty
        res = false;
        scaling = Inf;
        cert = true;
    end
    return
end
% cZ is not empty
        
% point or point cloud in constrained zonotope containment
if isnumeric(S)
    [res, cert, scaling] = aux_containsPoint(cZ,S,method,tol,scalingToggle);
    return
end

if isa(S, 'taylm')
    % Transform to pZ
    S = polyZonotope(S);
end

% Now we know that S is a contSet; again, the simpler this set is, the
% better the next few algorithms might run
S = compact(S);

% First, deal with trivial cases
if ~isBounded(S)
    % Unbounded -> not contained, since cZ is always
    % bounded
    res = false;
    cert = true;
    scaling = Inf;
    return
elseif representsa(S, 'emptySet')
    % Empty -> always contained
    res = true;
    cert = true;
    scaling = 0;
    return
else
    try
        [isPoint,p] = representsa(S, 'point');
        if isPoint
            [res, cert, scaling] = aux_containsPoint(cZ,p,method,tol,scalingToggle);
            return
        end
    catch ME
        if strcmp('CORA:notSupported', ME.identifier) || strcmp('MATLAB:maxlhs',ME.identifier)
            % If the code above returns an error either because there are
            % too many outputs, or the operation is not supported, we
            % conclude that it is not implemented, and we do nothing
        else
            % In any other case, something went wrong. Relay that
            % information.
            rethrow(ME);
        end
    end
end

% Now that the trivial cases have been dealt with, we need to choose
% which algorithm to use.

switch method
    case {'exact', 'exact:venum', 'exact:polymax'} % Exact algorithms
        [res, cert, scaling] = aux_exactParser(cZ, S, method, tol, maxEval, certToggle, scalingToggle);
    case {'approx', 'approx:st'} % Approximative algorithms
        [res, cert, scaling] = aux_approxParser(cZ, S, method, tol, maxEval, certToggle, scalingToggle);
    case {'sampling', 'sampling:primal', 'sampling:dual'} % Stochastic algorithms
        if ~isa(S, 'conZonotope')
            % For now, we only support the cases where S is a conZonotope
            throw(CORAerror('CORA:noSpecificAlg',method,cZ,S));
        end
        switch method
            case {'sampling', 'sampling:primal'}
                [res, cert, scaling] = aux_sampling(cZ, S, tol, maxEval, certToggle, scalingToggle);
            case 'sampling:dual'
                [res, cert, scaling] = aux_samplingDual(cZ, S, tol, maxEval, certToggle, scalingToggle);
        end
    otherwise
        throw(CORAerror('CORA:noSpecificAlg',method,cZ,S));
end
end


% Auxiliary functions -----------------------------------------------------

function [res, cert, scaling] = aux_exactParser(cZ, S, method, tol, maxEval, certToggle, scalingToggle)
    % Chooses what algorithm to call upon to check for exact containment
    switch class(S)
        % If S is convex, there are two main options: Either S is a
        % polyhedron, or it is not.
        % In the first case, it makes sense to compute the vertices
        % of S, in the latter we have to compute the halfspaces of
        % cZ.
        % The reason why in the first case we don't compute the
        % halfspace is simple: As of now, the
        % polytope-representation of cZ is computed by first
        % computing the vertices of cZ (expensive), and then deduce
        % the halfspace representation (even more expensive).
        case {'conHyperplane', 'emptySet', 'fullspace', ...
                'halfspace', 'interval', 'polytope', ...
                'zonoBundle', 'conZonotope', 'zonotope'}
            % Choose venum, unless the user decides otherwise
            if strcmp(method, 'exact')
                method = 'exact:venum';
            end

        case {'capsule', 'ellipsoid'}
            % Choose polymax, unless the user decides otherwise
            if strcmp(method, 'exact')
                method = 'exact:polymax';
            elseif strcmp(method, 'exact:venum')
                % Actually nevermind, venum is not an option here
                throw(CORAerror('CORA:noSpecificAlg',method,cZ,S));
            end
        otherwise
            throw(CORAerror('CORA:noExactAlg',cZ,S));
    end

    % Now, it only remains to execute the corresponding algorithm
    switch method
        case 'exact:venum'
            V = vertices(S);
            [res_V, ~, scaling_V] = contains_(cZ,V,'exact',tol,maxEval,certToggle,scalingToggle);
            res = all(res_V);
            cert = true;
            scaling = max(scaling_V);
        case 'exact:polymax'
            P_cZ = polytope(cZ);
            [res, cert, scaling] = contains_(P_cZ,S,'exact',tol,maxEval,certToggle,scalingToggle);
        otherwise
            throw(CORAerror('CORA:noSpecificAlg',method,cZ,S));
    end

end

function [res, cert, scaling] = aux_approxParser(cZ, S, method, tol, maxEval, certToggle, scalingToggle)
    % Chooses what algorithm to call upon to check for approx containment

    switch class(S)
        % If S is convex, there are two main options: Either S is a
        % polyhedron, or it is not.
        % In the first case, we can use the algorithm from Theorem 1 in
        % [1], otherwise we first need to transform the set into a
        % polyhedron (right now, this is done by computing the interval
        % over-approximation)
        case {'conHyperplane', 'emptySet', 'fullspace', ...
                'halfspace', 'interval', 'polytope', ...
                'zonoBundle', 'conZonotope', 'zonotope'}
            % Those are all the cases that can represent polyhedra
            S = conZonotope(S);
            res = aux_containsAHpolytope_approx(cZ,S,tol);
            cert = res;
            scaling = NaN;
            if scalingToggle
                throw(CORAerror('CORA:notSupported',...
                    "The computation of the scaling factor for the 'approx' method is not yet implemented for "+class(S)+"-in-"+class(cZ)+" containment."));
            end
        case {'capsule', 'ellipsoid'}
            % All other convex set representations
            if ~strcmp(method, 'approx')
                % This is to avoid people using st here, as the algorithm
                % that is being used here is not quite st.
                throw(CORAerror('CORA:noSpecificAlg',method,cZ,S));
            end
            
            cZ_S = conZonotope(interval(S));

            res = aux_containsAHpolytope_approx(cZ,cZ_S,tol);
            cert = res;
            scaling = NaN;
            if scalingToggle
                throw(CORAerror('CORA:notSupported',...
                    "The computation of the scaling factor for the 'approx' method is not yet implemented for "+class(S)+"-in-"+class(cZ)+" containment."));
            end
           
        otherwise
            % If S is not a convex set representation, our best chance is
            % to check for containment using the halfspace representation
            % of cZ

            if ~strcmp(method, 'approx')
                % This is to avoid people using st here, as the algorithm
                % that is being used here is not st at all.
                throw(CORAerror('CORA:noSpecificAlg',method,cZ,S));
            end
            
            if isa(S, 'taylm') || isa(S, 'polyZonotope')
                % For those set representations, the support function is
                % implemented
                P_cZ = polytope(cZ);
                [res, cert, scaling] = contains_(P_cZ,S,'approx',tol,maxEval,certToggle,scalingToggle);
            else
                try
                    S = interval(S);
                catch ME
                    throw(CORAerror('CORA:noSpecificAlg','approx',cZ,S));
                end
                [res, cert, scaling] = cZ.contains_(S,'exact',tol,maxEval,certToggle,scalingToggle);
                cert = res;
            end
    end   
end

function [res, cert, scaling] = aux_containsPoint(cZ,p, method,tol,scalingToggle)

    if strcmp(method, 'exact:polymax')
        % If the user insists on polymax, use that method
        P = polytope(cZ);
        [res,cert,scaling] = P.contains_(p, method, tol, 0, false, scalingToggle);
        return
    end

    % If the algorithm is not polymax, we can use venum
    cert = true(1,size(p,2));

    % use linear programming to check if a point is located inside a
    % constrained zonotope; see (20) in [2]


    % get object properties
    c = cZ.c;
    G = cZ.G;
    A = cZ.A;
    b = cZ.b;

    k = size(A,1);

    N = size(p,2);

    if ~scalingToggle
        % Checking whether p \in conZonotope(c,G,A,b) is equivalent to
        % checking whether [p;b] \in zonotope([G;A],[c;zeros([k 1])])
        p_ext = [p;kron(b,ones([1 N]))];
        Z_ext = zonotope([c;zeros([k 1])], [G;A]);
        [res, cert, scaling] = Z_ext.contains_(p_ext, method, tol, 0, false, false);
        return
    else
        cZ_center = center(cZ);
        % We are searching for the minimal delta, such that
        % p \in delta * (cZ - cZ_center) + cZ_center
        % Dividing by delta, and rearranging terms, we get the equivalent
        % problem of finding the maximal lambda (with lambda = 1/delta)
        % such that
        % lambda * (p-cZ_center) + cZ_center \in cZ

        n = size(G, 1);
        m = size(G, 2);

        % construct inequality constraints
        Aineq = [eye(m) zeros([m 1]);-eye(m) zeros([m 1]); zeros([1 m]) -1];
        bineq = [ones(2*m,1);0];
        
        % construct equality constraints
        %Aeq = [A zeros([k 1]);G -(p-cZ_center)]; % Just for reference;
        % will be executed later
        beq = [b;cZ_center-c];
        
        
        % solve linear program
        f = [zeros(m,1);-1];

        res = false([1 N]);
        scaling = Inf([1 N]);

        % init linprog struct
        problem.f = f;
        problem.Aineq = Aineq;
        problem.bineq = bineq;
        problem.beq = beq;
        problem.lb = [];
        problem.ub = [];

        for i=1:N
            Aeq = [A zeros([k 1]);G -(p(:,i)-cZ_center)];
            problem.Aeq = Aeq;

            [~,fval,exitflag] = CORAlinprog(problem);
        
            if exitflag == -2
                % LP is not feasible -> cZ must be degenerate, and p is not in
                % the same subspace
                res(i) = false;
                scaling(i) = Inf;
            elseif exitflag == -3
                % LP is unbounded -> p must be the center of cZ
                res(i) = true;
                scaling(i) = 0;
            elseif exitflag == 1
                scaling(i) = abs(1/fval);
                res(i) = scaling(i) <= 1+tol;
            else
                throw(CORAerror('CORA:solverIssue', 'linprog'));
            end
        end
    end

end

function res = aux_containsAHpolytope_approx(cZ1,cZ2, tol)
% check polytope in polytope containment according to Theorem 1 in [1]

    % convert to AH polytopes
    [Y,y,Hy,hy] = priv_AHpolytope(cZ1);
    [X,x,Hx,hx] = priv_AHpolytope(cZ2);
    
    nx = size(X,2);
    ny = size(Y,2);
    
    qx = length(hx);
    qy = length(hy);
    
    d = length(x);
    
    % construct constraint X = Y * T
    temp = repmat({Y},[1,nx]);
    A1 = [blkdiag(temp{:}),zeros(d*nx,qy*qx),zeros(d*nx,ny)];
    b1 = reshape(X,[d*nx,1]);

    % construct constraint y-x = Y * beta
    A2 = [zeros(d,nx*ny),zeros(d,qx*qy),Y];
    b2 = y-x;

    % construct constraint lambda * Hx = Hy * T
    Atemp = [];
    for j = 1:nx
        h = Hx(:,j);
        temp = repmat({h'},[1,qy]);
        Atemp = [Atemp;blkdiag(temp{:})];
    end

    temp = repmat({Hy},[1,nx]);
    A3 = [blkdiag(temp{:}),-Atemp,zeros(size(Atemp,1),ny)];
    b3 = zeros(size(A3,1),1);

    % construct overall equality constraints
    Aeq = [A1;A2;A3];
    beq = [b1;b2;b3];

    % construct constraint lambda * hx <= hy + Hy beta
    temp = repmat({hx'},[1,qy]);
    A1 = [zeros(qy,ny*nx),blkdiag(temp{:}),-Hy];
    b1 = hy;

    % construct constraint lambda >= 0
    A2 = [zeros(qx*qy,nx*ny),-eye(qx*qy),zeros(qx*qy,ny)];
    b2 = zeros(qx*qy,1);
    
    % constuct overall inequality constraints
    A = [A1;A2];
    b = [b1;b2];
    
    % add slack variables
    A = [A, [-eye(qy);zeros(qy*qx,qy)]; zeros(qy,size(A,2)) -eye(qy)];
    b = [b;zeros(qy,1)];
    
    Aeq = [Aeq,zeros(size(Aeq,1),qy)];
    
    % solve linear program
    f = [zeros(nx*ny+qy*qx+ny,1);ones(qy,1)];

    problem.f = f';
    problem.Aineq = A;
    problem.bineq = b;
    problem.Aeq = Aeq;
    problem.beq = beq;
    problem.lb = [];
    problem.ub = [];
    
    [val,~,exitflag] = CORAlinprog(problem); 
    
    % check for containment
    res = true;
    
    if exitflag < 0 || any(val(nx*ny+qy*qx+ny+1:end) > tol)
        res = false;
    end
end

function [res, cert, scaling] = aux_sampling(cZ, S, tol, maxEval, certToggle, scalingToggle)
    % Solves the containment problem using a similar technique to [3] for
    % zonotopes (using the vertex sampling algorithm)

    % Generate points uniformly on S
    p = S.randPoint_(maxEval,'uniform');
    
    % Compute center; we will then draw a ray from the center to the
    % border, passing through the points in p, in order to generate points
    % on the boundary.

    c = center(S);
    p = aux_boundaryPoint(S, p, c); % For the moment, this works only if S
    % is a constrained zonotope. In the future, this could be extended to
    % other convex set representations


    [res, cert, scaling] = cZ.contains_(p,'exact',tol,0,certToggle,scalingToggle);

    res = all(res); 
    if res
        cert = false;
    else
        cert = true;
    end
    scaling = max(scaling);

end

function [res, cert, scaling] = aux_samplingDual(cZ, S, tol, maxEval, certToggle, scalingToggle)
    % Solves the containment problem using a similar technique to [3] for
    % zonotopes (using the halfspace sampling algorithm)

    % Center the cZ and s so that cZ is somewhat centered at the origin
    c = center(cZ);
    cZ = cZ - c;
    S = S - c;
    

    % Generate points uniformly on cZ
    p = cZ.randPoint_(maxEval,'uniform');
    
    % Compute center; we will then draw a ray from the center to the
    % border, passing through the points in p, in order to generate points
    % on the boundary.

    c = center(cZ);
    p = aux_boundaryPoint(cZ, p, c);

    % get object properties
    c = cZ.c;
    G = cZ.G;
    cZ_A = cZ.A;
    cZ_b = cZ.b;

    n = size(G, 1);
    m = size(G, 2);
    k = size(cZ_A, 1);

    % We need to check whether the cZ is degenerate; the algorithm
    % changes a wee little bit if that is the case.
    [cZ_isFullDim, X] = isFullDim(cZ);
    
    if ~cZ_isFullDim
        X = [X zeros([n n-size(X,2)])];
    end

    res = true;
    cert = false;
    scaling = 0;

    for i = 1:maxEval
        q = p(:,i);

        % define parameters for linear program, to compute the normal
        % vector
        if cZ_isFullDim
            f = [zeros([n+k+2*m 1]);-q];
            A = [zeros([2*m n+k]) -eye(2*m) zeros([2*m n])];
            A = [A; c' cZ_b' ones([1 m]) ones([1 m]) zeros([1 n])];
            b = [zeros([2*m 1]); 1];
            Aeq = [eye(n) zeros([n k]) zeros([n m]) zeros([n m]) eye(n)];
            Aeq = [Aeq; -G' cZ_A' eye(m) -eye(m) zeros([m n])];
            beq = [zeros([n+m 1])];

            % init linprog struct
            problem.f = f;
            problem.Aineq = A;
            problem.bineq = b;
            problem.Aeq = Aeq;
            problem.beq = beq;
            problem.lb = [];
            problem.ub = [];

            % get unit normal vector of hit facet
            linOut = CORAlinprog(problem);
            s = linOut(n+k+2*m + 1:end);
        else
            % If cZ is degenerate, we have to ensure that the normal
            % vector is within the subspace in which cZ lives
            f = [zeros([n 1]);zeros([n+k+2*m 1]);-q];
            A = [zeros([2*m n]) zeros([2*m n+k]) -eye(2*m) zeros([2*m n])];
            A = [A; zeros([1 n]) c' cZ_b' ones([1 m]) ones([1 m]) zeros([1 n])];
            b = [zeros([2*m 1]); 1];
            Aeq = [zeros(n) eye(n) zeros([n k]) zeros([n m]) zeros([n m]) eye(n)];
            Aeq = [Aeq; zeros([m n]) -G' cZ_A' eye(m) -eye(m) zeros([m n])];
            beq = [zeros([n+m 1])];

            Aeq = [Aeq;X zeros(n) zeros([n k]) zeros([n m]) zeros([n m]) -eye(n)];
            beq = [beq; zeros([n 1])];

            % init linprog struct
            problem.f = f;
            problem.Aineq = A;
            problem.bineq = b;
            problem.Aeq = Aeq;
            problem.beq = beq;
            problem.lb = [];
            problem.ub = [];

            % get unit normal vector of hit facet
            linOut = CORAlinprog(problem);
            s = linOut(n+n+k+2*m + 1:end);
        end

        s = s / norm(s);

        L = supportFunc(S,s,'upper');

        scaling = max([scaling L]);

        if scaling > 1+tol
            res = false;
            cert = true;
            if ~scalingToggle
                return
            end
        end

    end

end

function p = aux_boundaryPoint(cZ, p, c_start)
    % Generates a point on the boundary of cZ, starting from c and passing
    % through p (does that multiple times if p is an array of points)
    N = size(p, 2);

    % Quick approximation of the ratio of cZ; this is assuming that the p
    % are uniformly distributed on cZ.
    radius_approx = sqrt(max(sum((p-c_start).^2)));

    % get object properties
    c = cZ.c;
    G = cZ.G;
    A = cZ.A;
    b = cZ.b;

    n = size(G, 1);
    m = size(G, 2);
    k = size(A, 1);

    % We need to check whether the cZ is degenerate; the algorithm
    % changes a wee little bit if that is the case.
    [cZ_isFullDim, X] = isFullDim(cZ);
    
    if ~cZ_isFullDim
        X = [X zeros([n n-size(X,2)])];
    end

    % construct inequality constraints
    Aineq = [eye(m) zeros([m 1]);-eye(m) zeros([m 1]); zeros([1 m]) -1];
    bineq = [ones(2*m,1);0];

    % construct equality constraints
    % (just for reference here)
    %Aeq = [A zeros([k 1]);G -d];
    %beq = [b;p-c];


    % cost function
    f = [zeros(m,1);-1];

    for i=1:N
        direction = p(:,i)-c_start;
        if norm(direction) <= 1e-4 * radius_approx
            % The point p(:,i) is quite close to c. Choosing a random
            % direction is a bit safer
            direction = randn([dim(cZ) 1]);

            direction = X * direction;
        end

        direction = direction ./ norm(direction);

        Aeq = [A zeros([k 1]);G -direction];
        beq = [b;p(:,i)-c];

        % init linprog struct
        problem.f = f;
        problem.Aineq = Aineq;
        problem.bineq = bineq;
        problem.Aeq = Aeq;
        problem.beq = beq;
        problem.lb = [];
        problem.ub = [];

        [val,fval,exitflag] = CORAlinprog(problem);

        p(:,i) = p(:,i) + abs(fval) * direction;
    end
end

% ------------------------------ END OF CODE ------------------------------
