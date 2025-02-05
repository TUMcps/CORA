function res = test_removeRedundantIds()
% test_removeRedundantIds - unit test function for removal of redundant elements
% from the ID-vector
%
% Syntax:
%    res = test_removeRedundantIds()
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

% basic test case
id = [1, 2, 1, 3, 2];
E = [1, 2; 3, 4; 5, 6; 7, 8; 9, 10];

[E_unique,id_unique] = removeRedundantIds(E,id);

% ... expected result
exp_E_unique = [6, 8; 12, 14; 7, 8];
exp_id_unique = [1, 2, 3];

% ... check result
assert(all(E_unique == exp_E_unique,'all'));
assert(all(id_unique == exp_id_unique));

% no redundant ids
id = [1, 2, 3];
E = [1, 2; 3, 4; 5, 6];

[E_unique,id_unique] = removeRedundantIds(E,id);

% ... expected result
exp_E_unique = E;
exp_id_unique = id;

% ... check result
assert(all(E_unique == exp_E_unique,'all'));
assert(all(id_unique == exp_id_unique));

% all ids are redundant
id = [1, 1, 1, 1];
E = [1, 2; 3, 4; 5, 6; 7, 8];

[E_unique,id_unique] = removeRedundantIds(E,id);

% ... expected result
exp_E_unique = sum(E,1);
exp_id_unique = 1;

% ... check result
assert(all(E_unique == exp_E_unique,'all'));
assert(all(id_unique == exp_id_unique));

% empty exponent matrix
id = [1, 2, 3];
E = [];

[E_unique,id_unique] = removeRedundantIds(E,id);

% ... expected result
exp_E_unique = [];
exp_id_unique = [1, 2, 3];

% ... check result
assert(all(E_unique == exp_E_unique,'all'));
assert(all(id_unique == exp_id_unique));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
