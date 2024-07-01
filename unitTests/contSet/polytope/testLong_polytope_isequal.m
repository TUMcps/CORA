function res = testLong_polytope_isequal
% testLong_polytope_isequal - unit test function for set equality check
%
% Syntax:
%    res = testLong_polytope_isequal()
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
% See also: none

% Authors:       Mark Wetzlinger
% Written:       07-June-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% number of tests and tolerance
nrTests = 50;
tol = 1e-6;

for i=1:nrTests
    
    % random dimension
    n = randi(20);
    nrCon = n+1 + randi([1, 10*n]);

    % only inequality constraints
    P = polytope.generateRandom('Dimension',n,'NrConstraints',nrCon);

    % 2. has to be equal to itself
    if ~isequal(P,P,tol)
        throw(CORAerror('CORA:testFailed'));
    end
    
    % 2. re-order constraints -> same polytope
    randOrder = randperm(nrCon);
    P_reordered = polytope(P.A(randOrder,:),P.b(randOrder));

    if ~isequal(P,P_reordered,tol)
        throw(CORAerror('CORA:testFailed'));
    end

    % 3. add n redundant constraints
    A_redundant = zeros(n,n); b_redundant = zeros(n,1);
    for j=1:n
        % compute support function in a random direction
        randDir = randn(n,1); randDir = randDir ./ vecnorm(randDir);
        val = supportFunc_(P,randDir,'upper');
        % generate redundant constraint by using a larger offset
        A_redundant(j,:) = randDir;
        b_redundant(j) = val + 0.01 + rand(1);
    end
    randOrder = randperm(nrCon);
    P_redundant = polytope([P.A(randOrder,:); A_redundant],...
        [P.b(randOrder); b_redundant]);
    
    if ~isequal(P,P_redundant,tol)
        throw(CORAerror('CORA:testFailed'));
    end

    % 4. add one irredundant constraint
    % compute support function in a random direction
    randDir = randn(n,1); randDir = randDir ./ vecnorm(randDir);
    val = supportFunc_(P,randDir,'upper');
    randOrder = randperm(nrCon);
    P_nonredundant = polytope([P.A(randOrder,:); randDir'],...
        [P.b(randOrder); val - 0.01 - rand(1)]);
    
    if isequal(P,P_nonredundant,tol)
        throw(CORAerror('CORA:testFailed'));
    end

end

% ------------------------------ END OF CODE ------------------------------
