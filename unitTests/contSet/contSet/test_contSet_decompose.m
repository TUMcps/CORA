function res = test_contSet_decompose
% test_contSet_decompose - unit test function of decomposition of a set
%    into several lower-dimensional projections
%
% Syntax:
%    res = test_contSet_decompose
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
% Written:       16-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% tolerance
tol = 1e-14;

% init several sets and partition
S = {interval([-2;1;4;0;2],[-2;5;6;1;3]); ...
     zonotope(zeros(5,1),eye(5)); ...
     ellipsoid(eye(5),ones(5,1)); ...
     polytope([eye(5); -ones(1,5)],ones(6,1)); ...
     fullspace(5); ...
     emptySet(5); ...
     capsule(zeros(5,1),ones(5,1),2); ...
    };
blocks = [1 2; 3 5];

% loop over all cases
for i=1:numel(S)

    % decompose
    S_dec = decompose(S{i},blocks);

    % true decomposition
    S_true = {project(S{i},blocks(1,1):blocks(1,2));...
              project(S{i},blocks(2,1):blocks(2,2))};

    assertLoop(all(arrayfun(@(j) isequal(S_dec{j},S_true{j},tol),...
                            1:2,'UniformOutput',true)),i);

end

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
