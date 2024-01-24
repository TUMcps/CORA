function res = test_levelSet_representsa
% test_levelSet_representsa - unit test function of representsa
%
% Syntax:
%    res = test_levelSet_representsa
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

% Authors:       Tobias Ladner
% Written:       17-January-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

% init levelSets ----------------------------------------------------------

% single equation
syms x y
eqs = x^2 + y^2 - 4;
ls1 = levelSet(eqs,[x;y],'==');

% multiple equations
eq1 = x^2 + y^2 - 4;
eq2 = x + y;
ls2 = levelSet([eq1;eq2],[x;y],{'<=';'<='});

% emptySet
ls_empty = levelSet.empty(2);

% fullspace
ls_Inf = levelSet.Inf(2);

% representsa emtpySet ----------------------------------------------------

resvec(end+1) = ~representsa(ls1,'emptySet');
resvec(end+1) = ~representsa(ls2,'emptySet');
resvec(end+1) = representsa(ls_empty,'emptySet');
resvec(end+1) = ~representsa(ls_Inf,'emptySet');


% representsa fullspace ---------------------------------------------------

resvec(end+1) = ~representsa(ls1,'fullspace');
resvec(end+1) = ~representsa(ls2,'fullspace');
resvec(end+1) = ~representsa(ls_empty,'fullspace');
resvec(end+1) = representsa(ls_Inf,'fullspace');


% gather results ----------------------------------------------------------
res = all(resvec);


% ------------------------------ END OF CODE ------------------------------
