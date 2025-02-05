function res = testLong_conZonotope_randPoint
% testLong_conZonotope_randPoint - unit test function of randPoint
%
% Syntax:
%    res = testLong_conZonotope_randPoint
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger, Adrian Kulmburg
% Written:       14-March-2021
% Last update:   22-May-2023 (AK, added all new methods)
%                04-February-2024 (AK, added non-degenerate cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
tol = 1e-9;

% check empty conZonotope object -> error
conZono = conZonotope.empty(1);
try 
    p = randPoint(conZono);
catch ME
    if ~strcmp(ME.identifier,'CORA:emptySet')
        rethrow(ME)
    end
end


% number of tests
nrOfTests = 100;

methods = {'standard', 'extreme', 'uniform', 'uniform:hitAndRun', 'uniform:ballWalk', 'uniform:billiardWalk'};

for i = 1:size(methods, 2)
    disp(methods{i})

for j=1:nrOfTests
    % random dimension
    n = randi([2,8]); % small because of containment checks
    
    % random center
    c = randn(n,1);
    % random number of generators
    m = n+randi(5);
    % random generator matrix
    G = randn(n,m);
    % make sure it's full rank
    while rank(G) < n
        G = randn(n,m);
    end

    % Make it a 4-cases analysis
    if mod(j,4) == 0
        % Case 1: cZ is a non-degenerate zonotope
        A = [];
        b = [];
    elseif mod(j,4) == 1
        % Case 2: cZ is a non-degenerate constrained zonotope
        P = polytope.generateRandom('Dimension', n);
        conZono = conZonotope(P);
        G = conZono.G;
        c = conZono.c;
        A = conZono.A;
        b = conZono.b;
    elseif mod(j,4) == 2
        % Case 3: random constraints so that conZonotope represents just a
        % point as A being diagional forces each independent factor to one
        % value
        A = diag(1+rand(m,1));
        b = sign(randn(m,1));
    elseif mod(j,4) == 3
        % Case 4: cZ is a degenerate zonotope
        % Create a constraint matrix that is surjective, to make sure it
        % has a solution
        k_constraints = randi([1,m-1]);
        A = randn([k_constraints m]);
        while rank(A) < k_constraints
            A = randn([k_constraints m]);
        end
        b = randn([k_constraints 1]);

        % random number of 'degeneracy'
        degeneracy = randi([1,n]);
        [U, Sigma, V] = svd(G);
        for j=1:degeneracy
            Sigma(n-j+1,n-j+1) = 0;
        end
        G = U * Sigma * V';
    end
    
    % instantiate conZonotope
    conZono = conZonotope(c,G, A, b);

    nrGens = size(G,2);
    
    % compute random points
    nrPts = 10;
    r = rng().Seed;
    rng(r * j); % Makes it easier to debug, keep this!
    p = randPoint(conZono,nrPts,methods{i});
    
    % check for containment in zonotope
    Z = zonotope(c,G);
    res = all(contains(Z,p,'exact',tol));
    
    % random constraints so that conZonotope represents just a point
    % as A being diagional forces each independent factor to one value
    A = diag(1+rand(nrGens,1));
    b = sign(randn(nrGens,1));
    % instantiate conZonotope with constraints
    conZono = conZonotope(c,G,A,b);
    
    % compute single possible point -> bug in randPoint
%    p = randPoint(conZono,1,methods{i});
    
    % compute point analytically
    p_analytical = zeros(n,1);
    conZonoG = conZono.G;
    beta = zeros(nrGens,1);
    for k=1:nrGens
        % constraints yield only one possible value for beta_k
        beta(k) = conZono.b(k) / conZono.A(k,k);
        
        % check if beta_k in [-1,1]
        assertLoop(abs(beta(k)) <= 1,j,k)
        
        % contribution of k-th generator
        kthGen = beta(k) * conZonoG(:,k);
        % add to point
        p_analytical = p_analytical + kthGen;
    end
    p = p_analytical;
    
    % compare results
    assertLoop(compareMatrices(p,p_analytical),i)
    
end
end

% ------------------------------ END OF CODE ------------------------------
