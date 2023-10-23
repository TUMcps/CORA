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
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger
% Written:       29-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% number of tests
nrTests = 25;

for i=1:nrTests
    
    % random dimension
    n = randi(5);

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

    % emptiness check
    if ~isBounded(P)
        res = false; return
    end

end

% ------------------------------ END OF CODE ------------------------------
