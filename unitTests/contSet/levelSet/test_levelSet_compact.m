function res = test_levelSet_compact
% test_levelSet_compact - unit test function of compact
%
% Syntax:
%    res = test_levelSet_compact
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
% Written:       25-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init symbolic variables
syms x y

% equations
eq1 = -x^2 - y^2 + 5;
eq2 = x^2 + y^2 - 5;

% if only one level set given, do nothing
ls = levelSet(eq1,[x;y],'<');
S = compact(ls);
res = isequal(ls,S);

% change of comparison operator
ls = levelSet([eq1,eq1],[x;y],{'<','<='});
S = compact(ls);
S_true = levelSet(eq1,[x;y],'<');
res(end+1,1) = isequal(S,S_true);

% level set composed of two inequalities
ls = levelSet([eq1;eq2],[x;y],{'<=','<='});

% compute reduced representation
S = compact(ls);

% true solution
S_true = levelSet(eq1,[x;y],'==');

% check
res(end+1,1) = isequal(S,S_true);


% level set composed of infeasible pair of inequalities
ls = levelSet([eq1;eq2],[x;y],{'<=','<'});

% compute reduced representation
S = compact(ls);

% true solution
S_true = emptySet(2);

% check
res(end+1,1) = isequal(S,S_true);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
