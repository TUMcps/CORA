function res = testLong_polytope_compact
% testLong_polytope_compact - unit test function of compact
%
% Syntax:
%    res = testLong_polytope_compact
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

% Authors:       Mark Wetzlinger
% Written:       05-December-2022
% Last update:   ---
% Last revision: 31-July-2023 (MW, rename '...compact')

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% number of tests
nrTests = 50;

% special method for 2D case
for i=1:nrTests

    nrCon = 4 + randi(26);

    % ensure polytope to be bounded
    A = [rand(1,2); -rand rand; -rand(1,2); -rand rand; randn(nrCon-4,2)];
    b = rand(nrCon,1);
    P = polytope(A,b);

    % compute minimal representation
    P_ = compact(P);

    % check random support function directions
    for j=1:5
        dir = randn(2,1);
        dir = dir ./ vecnorm(dir);
        % evaluate support function
        val1 = supportFunc(P_,dir);
        val2 = supportFunc(P_,dir);
        % compare values (second check for -/+Inf)
        if ~withinTol(val1,val2,1e-10) && ~(val1 == val2)
            res = false; return
        end
    end

end

for i=1:nrTests

    % random dimension
    n = randi(10);

    % init interval
    I = interval.generateRandom('Dimension',n);

    % convert to polytope (already in minimal H-representation)
    P = polytope(I);

    % rotate (does not create redundancies)
    [M,~,~] = svd(randn(n));
    P = M*P;

    % compute minimal H-representation
    P_ = compact(P);

    % both should have same number of constraints
    if length(P.b) ~= length(P_.b)
        res = false; return
    end


    % init box outside of smaller box
    A = [eye(n);-eye(n);eye(n);-eye(n)];
    b = [ones(2*n,1);2*ones(2*n,1)];
    P = polytope(A,b);

    % rotate (does not create redundancies)
    [M,~,~] = svd(randn(n));
    P = M*P;
    P = normalizeConstraints(P,'A');

    % compute minimal H-representation
    P_ = compact(P);

    % minimal H-representation should be first half of the redundant
    % representation
    if length(P_.b) ~= round(0.5*length(P.b)) ...
            || ~compareMatrices([P.A(1:2*n,:),P.b(1:2*n)]',[P_.A,P_.b]',1e-14)
        res = false; return
    end


    % init random polytope
    P = polytope.generateRandom('Dimension',n);
    
    % construct box around polytope
    B = box(P);
    % add those constraints to the original polytope (all redundant since
    % the original polytope is bounded, slight enlargement for stability)
    P_box = polytope([P.A;B.A],[P.b;B.b+1e-4]);

    % compute minimal H-representation
    P_ = compact(P_box);

    % should have at least 2*n less constraints than original polytope
    if length(P_.b) > length(P_box.b) - 2*n
        res = false; return
    end

end

% ------------------------------ END OF CODE ------------------------------
