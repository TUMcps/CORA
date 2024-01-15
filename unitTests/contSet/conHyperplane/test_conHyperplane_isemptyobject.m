function res = test_conHyperplane_isemptyobject
% test_conHyperplane_isemptyobject - unit test function of isemptyobject
%
% Syntax:
%    res = test_conHyperplane_isemptyobject
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
% Written:       03-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% empty set instantiation
hyp = conHyperplane.empty(2);
res(end+1,1) = ~isemptyobject(hyp);

% without constraints
a = [2 1]; b = 1;
hyp = conHyperplane(a,b);
res(end+1,1) = ~isemptyobject(hyp);

% with constraints
a = [3 2 -1]; b = 0.5;
C = [0 0 0]; d = 1;
hyp = conHyperplane(a,b,C,d);
res(end+1,1) = ~isemptyobject(hyp);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
