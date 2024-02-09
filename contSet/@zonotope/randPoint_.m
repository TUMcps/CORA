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
%           'uniform:hitAndRun', 'uniform:billiardWalk')
%
% Outputs:
%    p - random point in R^n
%
% Example: 
%    Z = zonotope([1;0],rand(2,5));
%    p = randPoint(Z);
% 
%    plot(Z); hold on;
%    scatter(p(1,:),p(2,:),16,'r');
%
% References:
%    [1] Robert L. Smith: Efficient Monte Carlo Procedures for Generating
%    Points Uniformly Distributed Over Bounded Regions, Operations Research
%    1984
%    [2] Boris T. Polyak, E. N. Gryazina: Billiard Walk - a New Sampling
%    Algorithm for Control and Optimization, IFAC, 2014
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/randPoint, interval/randPoint

% Authors:       Matthias Althoff, Mark Wetzlinger, Adrian Kulmburg, Severin Prenitzer
% Written:       23-September-2008 
% Last update:   25-June-2021 (MP, add type gaussian)
%                19-August-2022 (MW, integrate standardized pre-processing)
%                22-May-2023 (AK, implemented uniform sampling)
%                20-January-2023 (TL, added radius method)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% get object properties
c = Z.c; G = Z.G; n = dim(Z);
        
% no generators: zonotope is just a point
if isempty(G) || ~any(any(G))
    % replicate center N times
    p = repmat(c,1,N); return
end
        
% generate different types of random points
if strcmp(type,'standard')
    % degeneracy or all-zero generators do not have any effect here
    
    % take random values for factors
    factors = -1 + 2*rand(size(G,2),N);
    % sample points
    p = c + G * factors;
    
% sampling of extreme random points
elseif strcmp(type,'extreme')
    
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
        temp = randPoint(Z,N,type);
        p(ind,:) = temp;
        p = c + S*p;
        return;
    end
    
    % remove redundant generators
    Z = compact_(Z,'zeros',eps);
    Z = compact_(Z,'all',1e-3);
    
    % compute number of zonotope vertices
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
        ind = randperm(q);
        V = V(:,ind);
        p = V(:,1:N);
        
    else
        % compute vertices and additional points on the boundary
        V = vertices(Z);
        N_ = N - size(V,2);
        V_ = aux_getRandomBoundaryPoints(Z,N_);
        p = [V,V_];
    end

elseif strcmp(type,'uniform')
    if representsa_(Z,'parallelotope',eps)
        % If Z is a parallelotope, we have a significantly easier solution
        p = aux_randPointParallelotopeUniform(Z, N);
    else
        % Default algorithm for the uniform distribution is the billiard walk,
        % since it gives more accurate results, despite a slightly longer
        % computation time
        p = aux_randPointBilliard(Z, N);
    end
elseif strcmp(type,'uniform:billiardWalk')
    p = aux_randPointBilliard(Z, N);
elseif strcmp(type,'uniform:hitAndRun')
    p = aux_randPointHitAndRun(Z, N);

elseif strcmp(type,'radius')
    p = aux_randPointRadius(Z,N);

else
    throw(CORAerror('CORA:noSpecificAlg',type,Z));
    
end

end


% Auxiliary functions -----------------------------------------------------

function V = aux_getRandomVertices(Z,N)
% generate random vertices

    n = dim(Z); m = size(Z.G,2);
    V = zeros(m,N); cnt = 1; G = Z.G;
 
    % loop until the desired number of vertices is achieved
    while cnt <= N
        
        % generate random zonotope face
        temp = randperm(m);
        ind = temp(1:n-1);
        Q = G(:,ind);
        c = ndimCross(Q);
        v = sign(c'*G)';
        
        % generate random vertex on the zonotope face
        while true
           v_ = v;
           v_(ind) = sign((-1 + 2*rand(n-1,1)));
           if ~ismember(v_',V','rows')
              V(:,cnt) = v_;
              cnt = cnt + 1;
              break;
           end
        end
    end
    
    % compute vertices
    V = Z.c + G*V;
end

function V = aux_getRandomBoundaryPoints(Z,N)
% generate random points on the zonotope vertices

    G = Z.G;
    n = dim(Z); m = size(G,2);
    V = zeros(m,N);
 
    % loop until the desired number of vertices is achieved
    for i = 1:N
        
        % generate random zonotope face
        temp = randperm(m);
        ind = temp(1:n-1);
        Q = G(:,ind);
        c = ndimCross(Q);
        r = rand();
        if r > 0.5
           c = -c; 
        end
        V(:,i) = sign(c'*G);
        
        % generate random point on the zonotope face
        V(ind,i) = -1 + 2*rand(n-1,1);
    end
    
    % compute vertices
    V = Z.c + G*V;
end

function q = aux_numberZonoVertices(Z)
% compute the number of zonotope vertices

    n = dim(Z); m = size(Z.G,2);
    D = zeros(n,m);
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
    % The idea is pretty simple: The matrix G can be decomposed using a
    % SVD as G = USV'. U is just a pure rotation, so it doesn't impact the
    % uniformity of the distribution.
    % So we can 'act' as if G = SV; since we know, where S has zeros, we
    % can then choose the direction vector d in such a way, that
    % d \in Im(S) is always guaranteed.
    % It then just suffices to multiply our results by U at the very end.
    if ~isFullDim(Z)
        r = rank(Z.G);
        [U,S,V] = svd(Z.G);
        G = S*V';
        Z = zonotope(zeros([dim(Z) 1]), G);
    else
        r = dim(Z);
        U = eye(dim(Z));
    end
    
    G = Z.G;
    n = size(G,1);
    m = size(G,2);
    
    p0 = randPoint(Z);
    p = zeros(n,N);

    % create option to suppress command line output of linprog
    suppressPrint = optimoptions('linprog', 'Display', 'off');

    % sample N points
    for i = 1:N
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
            
        % define parameters for linear programs
        f = vertcat(1, zeros(m,1));
        A = horzcat(zeros(2*m,1), vertcat(eye(m), -eye(m)));
        b = ones(2*m,1);
        Aeq = horzcat(d, -G);
        beq =  - p0;

        % execute linear programs
        minZ = linprog(f,A,b,Aeq,beq,[],[],suppressPrint);
        maxZ = linprog(-f,A,b,Aeq,beq,[],[],suppressPrint);

        % sample line segment uniformly
        minC = minZ(1);
        maxC = maxZ(1);
        samplePoint = p0 + d * (minC + unifrnd(0,1) * (maxC - minC));
        % return sample point
        p(:,i) = samplePoint;
        p0 = samplePoint;
    end
    
    p = U*p + c;
end


function p = aux_randPointBilliard(Z,N) 

% m - number of generators, p0 - start point, tau - variable 
    %     influencing trajectory length, R0 - maximum reflection number, 
    %     p - return matrix
    c = Z.c;
    Z = Z - c;
    
    % We need to check whether the zonotope Z is degenerate; the algorithm
    % changes a wee little bit if that is the case.
    % The idea is pretty simple: The matrix G can be decomposed using a
    % SVD as G = USV'. U is just a pure rotation, so it doesn't impact the
    % uniformity of the distribution.
    % So we can 'act' as if G = SV; since we know, where S has zeros, we
    % can then choose the direction vector d in such a way, that
    % d \in Im(S) is always guaranteed.
    % It then just suffices to multiply our results by U at the very end.
    if ~isFullDim(Z)
        r = rank(Z.generators);
        [U,S,V] = svd(Z.generators);
        G = S*V';
        Z = zonotope(zeros([dim(Z) 1]), G);
    else
        r = dim(Z);
        U = eye(dim(Z));
        G = Z.generators;
    end
    
    
    %G = Z.generators;
    n = size(G,1);
    m = size(G,2);
    p0 = randPoint(Z);
    tau = norm(Z,2,'ub');
    R0 = 10 * dim(Z);
    
    p = zeros(n,N);

    % create option to suppress command line output of linprog
    suppressPrint = optimoptions('linprog', 'Display', 'off');
    
    
    % sample N points
    for i = 1:N
        % set segment start point, trajectory length and maximum remaining reflections
        q0 = p0;
        l = tau * -log(unifrnd(0,1));
        R = R0;

        % sample direction
        d = zeros(n,1);
        for j = 1:n
            if j > r
                % Cut out the parts of d that would be canceled by G (i.e., S*V')
                d(j) = 0;
            else
                d(j) = normrnd(0,1);
            end
        end
        

        while true
            % define parameters for linear program
            f = vertcat(1, zeros(m,1));
            A = horzcat(zeros(2*m,1), vertcat(eye(m), -eye(m)));
            b = ones(2*m,1);
            Aeq = horzcat(d, -G);
            beq = -q0;
            
            % get intersection point of direction vector with zonotope boundary
            maxZ = linprog(-f,A,b,Aeq,beq,[],[],suppressPrint);
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
                for j = 1:n
                    if j > r
                        % Cut out the parts of d that would be canceled by G (i.e., S*V')
                    d(j) = 0;
                    else
                        d(j) = normrnd(0,1);
                    end
                end
                continue;
            else
                R = R - 1;
            end
            
            % define parameters for linear program
            f = horzcat(zeros(1,m), -transpose(q));
            A = horzcat(vertcat(ones(1,m), -eye(m), -eye(m)), ...
                vertcat(zeros(1,n), transpose(G), -transpose(G)));
            b = vertcat(1, zeros(2*m,1));
        
            % get unit normal vector of hit facet
            linOut = linprog(f,A,b,[],[],[],[],suppressPrint);
            s = linOut(m+1:m+n);
            s = s / norm(s);
    
            % change start point, direction and remaining length
            q0 = q;
            d = d - 2 * dot(d,s) * s;
            if r ~= dim(Z)
                d(r+1:end) = zeros([dim(Z) - r 1]);
            end
            l = l - segmentLength;
        end
    end
    p = U * p + c;
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
