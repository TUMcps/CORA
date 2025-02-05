function p = randPoint_(cZ,N,type,varargin)
% randPoint_ - generates a random point inside a constrained zonotope
%
% Syntax:
%    p = randPoint_(cZ)
%    p = randPoint_(cZ,N)
%    p = randPoint_(cZ,N,type)
%    p = randPoint_(cZ,'all','extreme')
%
% Inputs:
%    cZ - conZonotope object
%    N - number of random points
%    type - type of the random point ('standard', 'extreme', 'uniform' or
%           'uniform:hitAndRun', 'uniform:billiardWalk',
%           'uniform:ballWalk')
%
% Outputs:
%    p - random point in the constrained zonotope
%
% Example: 
%    Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
%    A = [1 1 1]; b = 1;
%    cZ = conZonotope(Z,A,b);
%    p = randPoint(cZ,100);
%
%    figure; hold on;
%    plot(cZ,[1,2],'r');
%    plot(p(1,:),p(2,:),'.k','MarkerSize',10);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% References:
%    [1] Robert L. Smith: Efficient Monte Carlo Procedures for Generating
%    Points Uniformly Distributed Over Bounded Regions, Operations Research
%    1984
%    [2] Boris T. Polyak, E. N. Gryazina: Billiard Walk - a New Sampling
%    Algorithm for Control and Optimization, IFAC, 2014
%
% See also: contSet/randPoint, zonotope/randPoint

% Authors:       Niklas Kochdumper, Severin Prenitzer, Adrian Kulmburg
% Written:       30-October-2020
% Last update:   19-August-2022 (MW, integrate standardized pre-processing)
%                27-March-2023 (MW, rename randPoint_)
%                22-May-2023 (SP, implementing uniform sampling)
%                04-February-2024 (AK, implemented ball walk)
% Last revision: --- 

% ------------------------------ BEGIN CODE -------------------------------

    % call zonotope method if no constraints are present
    if isempty(cZ.A)
        p = randPoint_(zonotope(cZ.c,cZ.G),N,type); return
    end

    [res_point, rep_point] = representsa(cZ, 'point');
    if res_point
        p = kron(rep_point, ones([1 N]));
        return
    end
    
    % return all extreme points 
    if ischar(N) && strcmp(N,'all')
        p = vertices(cZ); return
    end
    
    % generate random points
    p = zeros(dim(cZ),N);
    
    if strcmp(type,'standard')
        for i = 1:N
            p(:,i) = aux_randPointStandard(cZ);
        end
    elseif strcmp(type,'extreme')
        for i = 1:N
            p(:,i) = aux_randPointExtreme(cZ);
        end
    % The default uniform sampler is hitAndRun, as otherwise we need to
    % transform it into a polytope
    elseif strcmp(type,'uniform') || strcmp(type,'uniform:hitAndRun')
        p = aux_randPointHitAndRun(cZ,N);
    elseif strcmp(type, 'uniform:billiardWalk')
        p = aux_randPointBilliard(cZ,N);
    elseif strcmp(type, 'uniform:ballWalk')
        p = aux_randPointBallWalk(cZ,N);
    else
        throw(CORAerror('CORA:noSpecificAlg',type,cZ));
    end
end


% Auxiliary functions -----------------------------------------------------

function p = aux_randPointStandard(cZ)    
% generate random point within the constrained zonotope

    % construct inequality constraints for the unit cube
    n = size(cZ.G,2);
    A = [eye(n);-eye(n)];
    b = [ones(n,1);ones(n,1)];

    % calculate null space of the constraints
    Neq = null(cZ.A);

    % calculate a single point that satisfies the constraints
    x0 = pinv(cZ.A)*cZ.b;

    % transform the constraints to the null space
    A_ = A*Neq;
    b_ = b-A*x0;

    % instantiate polytope
    P = polytope(A_,b_);

    % compute Chebychev center in the zonotope-factor null-space
    p = randPoint_(P,1,'standard');

    % convert center back to the normal zonotope factor space
    p_ = Neq*p + x0;

    % compute center of the constraint zonotope using the factors from the
    % Chebychev center in the factor space
    p = cZ.c + cZ.G * p_;
end

function p = aux_randPointExtreme(cZ)
% generate random point on boundary of a constrained zonotope

    % center constrained zonotope at origin
    c = center(cZ);
    cZ = cZ + (-c);

    % select random direction
    n = dim(cZ);
    d = rand(n,1) - 0.5*ones(n,1);
    d = d./norm(d);

    % compute farest point in this direction that is still in set
    [~,x] = supportFunc_(cZ,d,'upper');
    p = x + c;
end

function p = aux_randPointHitAndRun(cZ,N)
    % m - number of generators, p0 - start point, p - return matrix
    
    c = cZ.c;
    
    cZ = cZ - c;
    
    % We need to check whether the constrained zonotope cZ is degenerate;
    % the algorithm changes a wee little bit if that is the case.
    [res, X] = isFullDim(cZ);

    n = dim(cZ);

    if ~res
        X = [X zeros([n n-size(X,2)])];
    end
    
    % G - generator matrix, n - dimension, c - center,
    % m - generators number, nc - linear equality constraints number,
    % p0 - start point, p - return matrix
    n = dim(cZ);
    c = cZ.c;
    G = cZ.G;
    m = size(G,2);
    nc = size(cZ.A,1);

    % Try to sample at least one point; if this leads to an error, this
    % will be caught by the global function randPoint.m (for example,
    % if cZ is empty)
    p0 = cZ.randPoint(1);

    p = zeros(n,N);

    % sample N points
    for i = 1:N
        % sample movement vector
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
        problem.Aeq = vertcat(horzcat(d, -G), horzcat(zeros(nc,1), cZ.A));
        problem.beq = vertcat(c - p0, cZ.b);
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

function p = aux_randPointBallWalk(cZ,N)
    n = dim(cZ);
    
    % First, make sure that the zonotope is zero-centered; save the center
    % for the end, when we need to translate back our points
    c = cZ.center;
    cZ = cZ - c;
    
    % We first try to 'round' the cZ as much as we can, by computing
    % the interval over-approximation
    I = interval(cZ);
    % Side lengths
    L = (I.sup - I.inf)./2;
    % If we have an unbounded polytope, ignore the unbounded parts
    L(abs(L) == Inf) = 1;
    % Same if we have a degenerate dimension
    L(abs(L) < 100 * eps) = 1;

    cZtrans = diag(1./L) * cZ;
    
    % Need to check whether the cZ is degenerate
    [res, X] = isFullDim(cZtrans);
    n = dim(cZtrans);
    if ~res
        X = [X zeros([n n-size(X,2)])];
    end

    % Radius of the ball
    delta = 1/sqrt(n);
    
    % Start with some random point, not necessarily uniformly distributed
    p0 = randPoint(cZtrans);
    p = zeros(n,N);
    % sample N points
    for i = 1:N
        while true
            %generate a random direction
            d = randn([dim(cZtrans) 1]);
            d = d./norm(d);
            
            % in case the polytope is degenerate, the direction can only be in
            % the subspace spanned by the vectors of X (which are normalized
            % and orthogonal, so all is well, we still have a uniformly sampled
            % direction)
            d = X * d;

            % Sample length
            ell = rand(1)^(1/n);

            % Build candidate point
            candidate = p0 +  delta * ell * d;
            
            if cZtrans.contains_(candidate,'exact',1e-5,0,false,false)
                p0 = candidate;
                p(:,i) = p0;
                break
            end
            
        end
    end
    
    % Transform everything back now
    p = diag(L) * p + c;
end


function p = aux_randPointBilliard(cZ,N)
    % Need to compute the center of cZ, to make sure that it contains the
    % origin
    center_cZ = center(cZ);
    cZ = cZ - center_cZ;

    n = dim(cZ);
    
    % We need to check whether the cZ is degenerate; the algorithm
    % changes a wee little bit if that is the case.
    [cZ_isFullDim, X] = isFullDim(cZ);
    
    if ~cZ_isFullDim
        X = [X zeros([n n-size(X,2)])];
    end
    
    
    c = cZ.c; % Careful, that's the center vector, NOT the center of cZ!
    G = cZ.generators;
    cZ_A = cZ.A;
    cZ_b = cZ.b;
    k = size(cZ_A,1);
    n = size(G,1);
    m = size(G,2);
    p0 = randPoint(cZ);
    tau = norm_(interval(cZ),2);
    R0 = 10 * dim(cZ);
    
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

        % in case the cZ is degenerate, the direction can only be in
        % the subspace spanned by the vectors of X (which are normalized
        % and orthogonal, so all is well, we still have a uniformly sampled
        % direction)
        d = X * d;
        

        while true
            % define parameters for linear program
            f = vertcat(1, zeros(m,1));
            A = horzcat(zeros(2*m,1), vertcat(eye(m), -eye(m)));
            b = ones(2*m,1);
            Aeq = [-d G; zeros([k 1]) cZ_A];
            beq =  [q0 - c; cZ_b];

            % init linprog struct
            problem.f = -f;
            problem.Aineq = A;
            problem.bineq = b;
            problem.Aeq = Aeq;
            problem.beq = beq;
            problem.lb = [];
            problem.ub = [];
            
            % get intersection point of direction vector with cZ boundary
            maxZ = CORAlinprog(problem);
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

            
            % sample new direction from last returned point if maximum 
            % reflection number is reached
            if R <= 0
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
    
            % change start point, direction and remaining length
            q0 = q;
            d = d - 2 * dot(d,s) * s;
            l = l - segmentLength;
        end
    end
    p = p + center_cZ;
end

% ------------------------------ END OF CODE ------------------------------
