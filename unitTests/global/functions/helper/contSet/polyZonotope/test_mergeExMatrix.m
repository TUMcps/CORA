function res = test_mergeExMatrix()
% test_mergeExMatrix - unit test function for merging of exponent matrices
%
% Syntax:
%    res = test_mergeExMatrix()
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

% Authors:       Lukas Schäfer
% Written:       07-February-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% identical ID-vectors
id1 = [1, 2, 3]';
id2 = [1, 2, 3]';
E1 = [1, 2; 3, 4; 5, 6];
E2 = [7, 8; 9, 10; 11, 12];

[id, E1_new, E2_new] = mergeExpMatrix(id1, id2, E1, E2);

% ... check results
assert(isequal(id, id1));
assert(isequal(E1_new, E1));
assert(isequal(E2_new, E2));

% partially overlapping ID-vectors
id1 = [1, 2, 3]';
id2 = [3, 4, 5]';
E1 = [1, 2; 3, 4; 5, 6];
E2 = [7, 8; 9, 10; 11, 12];

[id, E1_new, E2_new] = mergeExpMatrix(id1, id2, E1, E2);

% ... expected results
exp_id =  [1, 2, 3, 4, 5]';
exp_E1 = [1, 2; 3, 4; 5, 6; 0, 0; 0, 0];
exp_E2 = [0, 0; 0, 0; 7, 8; 9, 10; 11, 12];

% ... check results
assert(isequal(id, exp_id));
assert(isequal(E1_new, exp_E1));
assert(isequal(E2_new, exp_E2));

% no common ids
id1 = [1, 2, 3]';
id2 = [4, 5, 6]';
E1 = [1, 2; 3, 4; 5, 6];
E2 = [7, 8; 9, 10; 11, 12];

[id, E1_new, E2_new] = mergeExpMatrix(id1, id2, E1, E2);

% ... expected results
exp_id = [1, 2, 3, 4, 5, 6]';
exp_E1 = [1, 2; 3, 4; 5, 6; 0, 0; 0, 0; 0, 0];
exp_E2 = [0, 0; 0, 0; 0, 0; 7, 8; 9, 10; 11, 12];

% ... check results
assert(isequal(id, exp_id));
assert(isequal(E1_new, exp_E1));
assert(isequal(E2_new, exp_E2));

% one empty ID-vector
id1 = [1, 2, 3]';
id2 = [];
E1 = [1, 2; 3, 4; 5, 6];
E2 = [];
    
[id, E1_new, E2_new] = mergeExpMatrix(id1, id2, E1, E2);

% ... expected results
exp_id =  id1;
exp_E1 = E1;

% ... check results
assert(isequal(id, exp_id));
assert(isequal(E1_new, exp_E1));
assert(isempty(E2_new));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
