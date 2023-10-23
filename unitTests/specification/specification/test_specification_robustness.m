function res = test_specification_robustness
% test_specification_robustness - unit test for robustness
%
% Syntax:
%    res = test_specification_robustness
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
% Written:       30-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% halfspace and point
hs = halfspace([1 0],1);
p = [0;0];

% unsafe set
spec = specification(hs,'unsafeSet');
% compute robustness
val = robustness(spec,p);
res = withinTol(val,-1);

% safe set
spec = specification(hs,'safeSet');
% compute robustness
val = robustness(spec,p);
res(end+1,1) = withinTol(val,1);

% invariant
spec = specification(hs,'invariant');
% compute robustness
val = robustness(spec,p);
res(end+1,1) = withinTol(val,1);

% multiple points
p = [0 0; 1 -1; 4 1; -8 -3]';
val = robustness(spec,p);
res(end+1,1) = all(withinTol(val,[1 0 -3 9]));


% other sets
I = interval([-2;-1;0],[1;4;2]);
spec = specification(I,'safeSet');
p = [0;0;0];
val = robustness(spec,p);
res(end+1,1) = withinTol(val,0);

Z = zonotope([1;-1],[1 2; -1 1]);
spec = specification(Z,'unsafeSet');
p = [2;2];
val = robustness(spec,p);
res(end+1,1) = withinTol(val,1);


% multiple specifications
Z = zonotope([1;-1],[1 2; -1 1]);
I = interval([-2;-1],[1;4]);
spec = [specification(Z,'safeSet');specification(I,'safeSet')];
p = [2;2];
val = robustness(spec,p);
res(end+1,1) = withinTol(val,-1);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
