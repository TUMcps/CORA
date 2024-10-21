function res = testLong_polytope_isBounded
% testLong_polytope_isBounded - unit test function for boundedness
%
% Syntax:
%    res = testLong_polytope_isBounded()
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
% Written:       29-November-2022
% Last update:   04-June-2024 (MW, test for unboundedness)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% number of tests
nrTests = 50;

for i=1:nrTests
    
    % random dimension
    n = randi(10);

    % compute random vertices
    I = interval(-ones(n,1),ones(n,1));
    temp = randPoint(I,100);
    % take extreme point in each dimension
    V = zeros(n,2*n);
    for j=1:n
        [~,maxIdx] = max(temp(j,:));
        [~,minIdx] = min(temp(j,:));
        V(:,(2*j-1)) = temp(:,maxIdx);
        V(:,2*j) = temp(:,minIdx);
        temp(:,[maxIdx,minIdx]) = [];
    end

    % init polytope
    P = polytope(V);

    % boundedness check
    assertLoop(isBounded(P),i)

    % sample normal vectors only from one halfspace -> polytope unbounded
    nrCon = randi([n+1,2*n]);
    randDir = randn(n,1);
    randDir = randDir / vecnorm(randDir);
    A = zeros(0,n);
    while size(A,1) < nrCon
        randDir_j = randn(n,1);
        % ensure that randDir_j looks into same halfspace as randDir
        if randDir' * randDir_j < 0
            continue
        end
        A = [A; randDir_j'];
    end
    % select random offsets (don't affect boundedness)
    b = rand(nrCon,1);

    % init polytope
    P = polytope(A,b);

    % boundedness check
    assertLoop(~isBounded(P),i)

end

% ------------------------------ END OF CODE ------------------------------
