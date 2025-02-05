function p = randPoint_(Z,N,type,varargin)
% randPoint_ - generates random points within a zonotope
%
% Syntax:
%    p = randPoint_(Z)
%    p = randPoint_(Z,N)
%    p = randPoint_(Z,N,type)
%    p = randPoint_(Z,'all','extreme')
%
% Inputs:
%    Z - zonotope object
%    N - number of random points
%    type - type of the random point ('standard', 'extreme', 'uniform' or
%           'uniform:hitAndRun', 'uniform:billiardWalk',
%           'uniform:ballWalk')
%
% Outputs:
%    p - random point (cloud) in R^n
%
% Example: 
%    Z = zonotope([1;0],[1 0 1; -1 2 1]);
%    p = randPoint(Z);
% 
%    plot(Z); hold on;
%    scatter(p(1,:),p(2,:),16,'r');
%
% References:
%    [1] Robert L. Smith: Efficient Monte Carlo Procedures for Generating
%        Points Uniformly Distributed Over Bounded Regions, Operations
%        Research, 1984.
%    [2] Boris T. Polyak, E. N. Gryazina: Billiard Walk - a New Sampling
%        Algorithm for Control and Optimization, IFAC, 2014
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/randPoint, interval/randPoint_

% Authors:       Matthias Althoff, Mark Wetzlinger, Adrian Kulmburg, Severin Prenitzer
% Written:       23-September-2008 
% Last update:   25-June-2021 (MP, add type gaussian)
%                19-August-2022 (MW, integrate standardized pre-processing)
%                22-May-2023 (AK, implemented uniform sampling)
%                20-January-2024 (TL, added radius method)
%                03-March-2024 (TL, made boundary method accessible)
% Last revision: 05-October-2024 (MW, refactor)

% ------------------------------ BEGIN CODE -------------------------------
        
% zonotope is just a point -> replicate center N times
if representsa_(Z,'point',eps)
    p = repmat(Z.c,1,N);
    return
end
        
% generate different types of random points
switch type
    case 'standard'
        p = aux_randPoint_standard(Z,N);

    case 'extreme'
        % sampling of extreme random points
        p = aux_randPoint_extreme(Z,N);

    case {'uniform','uniform:billiardWalk'}
        if representsa_(Z,'parallelotope',eps)
            % for parallelotopes, there is a significantly easier solution
            p = aux_randPointParallelotopeUniform(Z,N);
        else
            % Default algorithm for the uniform distribution is the
            % billiard walk, since it gives more accurate results, despite
            % a slightly longer computation time
            p = aux_randPointBilliard(Z,N);
        end
    case 'uniform:ballWalk'
        p = aux_randPointBallWalk(Z, N);
    case 'uniform:hitAndRun'
        p = aux_randPointHitAndRun(Z, N);
    case 'radius'
        p = aux_randPointRadius(Z,N);
    case 'boundary'
        p = aux_getRandomBoundaryPoints(Z,N);
    otherwise
        throw(CORAerror('CORA:noSpecificAlg',type,Z));
    
end

end


% Auxiliary functions -----------------------------------------------------

function p = aux_randPoint_standard(Z,N)
% degeneracy or all-zero generators do not have any effect here
    
    % take random values for factors
    factors = -1 + 2*rand(size(Z.G,2),N);
    % sample points
    p = Z.c + Z.G * factors;

end

function p = aux_randPoint_extreme(Z,N)

n = dim(Z);
c = Z.c; G = Z.G;

% 1D case
if n == 1
    % flush all generators into one
    G = sum(abs(G),2);
    % random signs
    s = sign(randn(1,N));
    % instantiate points
    p = c + s*G;
    return
end

% remove redundant generators
Z = compact_(Z,'all',1e-10);

% consider degenerate case
if rank(G) < n
    Z = Z + (-c);
    p = c * zeros(1,N);
    [S,V,~] = svd([-G,G]);
    d = diag(V);
    ind = find(d > eps);
    if isempty(ind)
        return;
    end
    Z = project(S'*Z,ind);
    p_proj = randPoint(Z,N,'extreme');
    p(ind,:) = p_proj;
    p = c + S*p;
    return;
end

% compute approximate number of zonotope vertices
q = aux_numberZonoVertices(Z);

if ischar(N) && strcmp(N,'all')
    % return all extreme points
    p = vertices(Z);
    
elseif 10*N < q
    % generate random vertices
    p = aux_getRandomVertices(Z,N);
    
elseif N < q
    % select random vertices
    V = vertices(Z);
    ind = randperm(size(V,2));
    V = V(:,ind);
    p = V(:,1:N);
    
else
    % compute vertices and additional points on the boundary
    V = vertices(Z);
    N_ = N - size(V,2);
    V_ = aux_getRandomBoundaryPoints(Z,N_);
    p = [V,V_];
end

end

function V = aux_getRandomVertices(Z,N)
% generate random vertices

    n = dim(Z); nrGen = size(Z.G,2);
    V = zeros(nrGen,N); cnt = 1;
 
    % loop until the desired number of vertices is achieved
    while cnt <= N
        
        % generate random zonotope face
        randOrder = randperm(nrGen);
        ind = randOrder(1:n-1);
        Q = Z.G(:,ind);
        c = ndimCross(Q);
        v = sign(c'*Z.G)';
        
        % generate random vertex on the zonotope face
        while true
           v_ = v;
           v_(ind) = sign(-1 + 2*rand(n-1,1));
           if ~ismember(v_',V','rows')
              V(:,cnt) = v_;
              cnt = cnt + 1;
              break;
           end
        end
    end
    
    % compute vertices
    V = Z.c + Z.G*V;
end

function V = aux_getRandomBoundaryPoints(Z,N)
% generate random points on the zonotope vertices

    n = dim(Z); nrGen = size(Z.G,2);
    V = zeros(nrGen,N);
 
    % loop until the desired number of vertices is achieved
    for i = 1:N
        
        % generate random zonotope face
        randOrder = randperm(nrGen);
        ind = randOrder(1:n-1);
        Q = Z.G(:,ind);
        c = ndimCross(Q);
        if rand() > 0.5
           c = -c; 
        end
        V(:,i) = sign(c'*Z.G);
        
        % generate random point on the zonotope face
        V(ind,i) = -1 + 2*rand(n-1,1);
    end
    
    % compute vertices
    V = Z.c + Z.G*V;
end

function q = aux_numberZonoVertices(Z)
% compute the number of zonotope vertices

    n = dim(Z); nrGen = size(Z.G,2);

    if nrGen == 0
        % only center
        q = 1;
        return;
    end

    D = zeros(n,nrGen);
    D(1,:) = 2*ones(1,size(D,2));
    D(:,1) = 2*ones(size(D,1),1);

    for j = 2:size(D,1)
        for k = 2:size(D,2)
            D(j,k) = D(j,k-1) + D(j-1,k-1);
        end
    end
    
    q = D(end,end); 
end

function p = aux_randPointHitAndRun(Z,N)
    % m - number of generators, p0 - start point, p - return matrix
    
    c = Z.c;
    Z = Z - c;
    
    % We need to check whether the zonotope Z is degenerate; the algorithm
    % changes a wee little bit if that is the case.
    [res, X] = isFullDim(Z);

    n = dim(Z);

    if ~res
        X = [X zeros([n n-size(X,2)])];
    end
    
    G = Z.G;
    m = size(G,2);
    
    p0 = randPoint(Z);
    p = zeros(n,N);

    % sample N points
    for i = 1:N
        d = randn([n 1]);
        d = d./norm(d);

        % in case the zonotope is degenerate, the direction can only be in
        % the subspace spanned by the vectors of X (which are normalized
        % and orthogonal, so all is well, we still have a uniformly sampled
        % direction)
        d = X * d;
            
        % define parameters for linear programs
        f = vertcat(1, zeros(m,1));
        problem.Aineq = horzcat(zeros(2*m,1), vertcat(eye(m), -eye(m)));
        problem.bineq = ones(2*m,1);
        problem.Aeq = horzcat(d, -G);
        problem.beq = -p0;
        problem.lb = [];
        problem.ub = [];

        % execute linear programs
        problem.f = f;
        minZ = CORAlinprog(problem);
        problem.f = -f;
        maxZ = CORAlinprog(problem);

        % sample line segment uniformly
        minC = minZ(1);
        maxC = maxZ(1);
        samplePoint = p0 + d * (minC + unifrnd(0,1) * (maxC - minC));
        % return sample point
        p(:,i) = samplePoint;
        p0 = samplePoint;
    end
    
    p = p + c;
end

function p = aux_randPointBallWalk(Z,N)
    n = dim(Z);
    
    % First, make sure that the zonotope is zero-centered; save the center
    % for the end, when we need to translate back our points
    c = Z.center;
    Z = Z - c;
    
    % Then, for the ball-walk, we need to make sure that the inscribed
    % ellipsoid in Z is actually a sphere. For this, we need the SVD:
    [U, S, V] = svd(Z.generators);
    % U is a rotation, which is bijective; we can thus turn Z:
    Z = U' * Z;
    % We now have a zonotope that has its main axes aligned with the
    % canonical basis. We compute the rank of S now, to check if the
    % zonotope is degenerate.
    % If not, we 'invert' the scaling introduced by S. If it is degenerate,
    % we only do this for the relevant directions:
    r = rank(S);
    s = diag(S);
    s = 1./s(1:r);
    s = [s;ones([n-r 1])];
    S_inv = diag(s);
    
    Z = S_inv * Z;
    
    % We have now turned and scaled Z in such a way, that it contains the
    % unit ball, and is contained in the unit ball scaled by sqrt(n), where
    % n is the dimension. We choose the length of the ball-walk to be
    % 1/sqrt(n)
    % this could be improved in the future, but at the moment it yields
    % pretty good results overall
    delta = 1/sqrt(n);
    
    % Start with some random point, not necessarily uniformly distributed
    p0 = randPoint(Z);
    p = zeros(n,N);

    % sample N points
    for i = 1:N
        while true
            % sample movement vector
            d = zeros(n,1);
            for j = 1:n
                if j > r
                    % Cut out the parts of d that would be canceled by G (i.e., S*V')
                    d(j) = 0;
                else
                    d(j) = normrnd(0,1);
                end
            end

            % Sample length
            ell = rand(1)^(1/n);

            % Build candidate point
            candidate = p0 +  delta * ell * d;
            
            if Z.contains_(candidate,'exact',1e-5,0,false,false)
                p0 = candidate;
                p(:,i) = p0;
                break
            end
            
        end
    end
    
    % Transform everything back now
    p = U * diag(1./s) * p + c;
end


function p = aux_randPointBilliard(Z,N) 
    % m - number of generators, p0 - start point, tau - variable 
    %     influencing trajectory length, R0 - maximum reflection number, 
    %     p - return matrix

    c = Z.c;
    Z = Z - c;

    n = dim(Z);
    
    % We need to check whether the zonotope Z is degenerate; the algorithm
    % changes a wee little bit if that is the case.
    [Z_isFullDim, X] = isFullDim(Z);

    if ~Z_isFullDim
        X = [X zeros([n n-size(X,2)])];
    end
    
    
    G = Z.generators;
    m = size(G,2);
    p0 = randPoint(Z);
    tau = norm(Z,2,'ub');
    R0 = 10 * dim(Z);
    
    p = zeros(n,N);    
    
    % sample N points
    for i = 1:N
        % set segment start point, trajectory length and maximum remaining reflections
        q0 = p0;
        l = tau * -log(unifrnd(0,1));
        R = R0;

        % sample direction
        d = randn([n 1]);
        d = d./norm(d);

        % in case the polytope is degenerate, the direction can only be in
        % the subspace spanned by the vectors of X (which are normalized
        % and orthogonal, so all is well, we still have a uniformly sampled
        % direction)
        d = X * d;
        

        while true
            % define parameters for linear program
            problem_boundary.f = vertcat(1, zeros(m,1));
            problem_boundary.Aineq = horzcat(zeros(2*m,1), vertcat(eye(m), -eye(m)));
            problem_boundary.bineq = ones(2*m,1);
            problem_boundary.Aeq = horzcat(d, -G);
            problem_boundary.beq = -q0;
            problem_boundary.lb = [];
            problem_boundary.ub = [];
            
            % get intersection point of direction vector with zonotope boundary
            maxZ = CORAlinprog(problem_boundary);
            maxC = maxZ(1);
            q = q0 + maxC * d;

            % return point on line segment if the trajectory length is reached
            lineSegment = q - q0;
            segmentLength = norm(lineSegment);
            
            if l <= segmentLength
                p(:,i) = q0 + (l / segmentLength) * lineSegment;
                p0 = p(:,i);
                break;
            end

            % check if a nonsmooth boundary point is met
            extremeCounter = 0;
            for j = 2:m+1
                if abs(abs(maxZ(j)) - 1) < eps(1)
                    extremeCounter = extremeCounter + 1;
                end
            end
            

            % sample new direction from last returned point if maximum 
            % reflection number or nonsmooth boundary point is reached
            if R <= 0 || extremeCounter <= m - 2 || extremeCounter == m
                R = R0;
                q0 = p0;
                l = tau * -log(unifrnd(0,1));

                d = randn([n 1]);
                d = d./norm(d);
                d = X * d;

                continue;
            else
                R = R - 1;
            end
            
            % define parameters for linear program
            if Z_isFullDim
                problem_facet.f = [zeros([m 1]); -q];
                problem_facet.Aineq = [ones([1 m]) zeros([1 n]);...
                                        -eye(m) G';...
                                        -eye(m) -G'];
                problem_facet.bineq = [1; zeros(2*m,1)];
                problem_facet.Aeq = [];
                problem_facet.beq = [];
                problem_facet.lb = [];
                problem_facet.ub = [];

                % get unit normal vector of hit facet
                linOut = CORAlinprog(problem_facet);
                s = linOut(m+1:m+n);
            else
                % Need to ensure that the normal vector is in the same
                % subspace as Z
                problem_facet.f = [zeros([n+m 1]); -q];
                problem_facet.Aineq = [zeros([1 n]) ones([1 m]) zeros([1 n]);...
                                       zeros([m n]) -eye(m) G';...
                                       zeros([m n]) -eye(m) -G'];
                problem_facet.bineq = [1; zeros(2*m,1)];
                problem_facet.Aeq = [X zeros([n m]) -eye(n)];
                problem_facet.beq = zeros([n 1]);
                problem_facet.lb = [];
                problem_facet.ub = [];
            
                % get unit normal vector of hit facet
                linOut = CORAlinprog(problem_facet);
                s = linOut(n+m+1:n+m+n);
            end

            s = s / norm(s);
    
            % change start point, direction and remaining length
            q0 = q;
            d = d - 2 * dot(d,s) * s;
            l = l - segmentLength;
        end
    end
    p = p + c;
end

function p = aux_randPointRadius(Z,N)

    % sample random points on boundary to obtain radii
    radii = aux_getRandomBoundaryPoints(Z,N) - Z.c;
    
    % sample semi-uniformly within 'sphere'
    p = Z.c + nthroot(rand(1,N),dim(Z)).*radii;

end

function p = aux_randPointParallelotopeUniform(Z, N)
    Z = compact(Z);

    G = generators(Z);
    c = center(Z);

    % Uniformly sample unit hypercube
    p = 2.*rand([size(G,2) N]) - 1;

    p = G*p + c;
end


% ------------------------------ END OF CODE ------------------------------
