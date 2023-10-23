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
%    type - type of the random point ('standard', 'extreme', 'hit-and-run'
%       or 'billiard-walk')
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
% Last revision: --- 

% ------------------------------ BEGIN CODE -------------------------------

    % call zonotope method if no constraints are present
    if isempty(cZ.A)
        p = randPoint_(zonotope(cZ.c,cZ.G),N,type); return
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
        % G - generator matrix, n - dimension, c - center, 
        % m - generators number, nc - linear equality constraints number, 
        % p0 - start point, p - return matrix
        n = dim(cZ);
        c = cZ.c;
        G = cZ.G;s
        m = size(G,2);
        nc = size(cZ.A,1);
        
        % Try to sample at least one point; if this leads to an error, this
        % will be caught by the global function randPoint.m (for example,
        % if cZ is empty)
        p0 = cZ.randPoint(1);
        
        p = zeros(n,N);
        
        % create option to suppress command line output of linprog
        persistent suppressPrint
        if isempty(suppressPrint)
            suppressPrint = optimoptions('linprog', 'Display', 'off');
        end
        
        % sample N points
        for i = 1:N
            % sample movement vector
            d = zeros(n,1);
            for j = 1:n
                d(j) = normrnd(0,1);
            end
            
            % define parameters for linear programs
            f = vertcat(1, zeros(m,1));
            A = horzcat(zeros(2*m,1), vertcat(eye(m), -eye(m)));
            b = ones(2*m,1);
            Aeq = vertcat(horzcat(d, -G), horzcat(zeros(nc,1), cZ.A));
            beq = vertcat(c - p0, cZ.b);
            
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
    elseif strcmp(type, 'uniform:billiardWalk')
        p = randPoint(polytope(cZ),N,'uniform:billiardWalk');
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

% ------------------------------ END OF CODE ------------------------------
