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

% instantiate constrained hyperplanes
hyp1 = conHyperplane();

a = [3; 2; -1];
b = 0.5;
C = [0 0 0];
d = 1;
hyp2 = conHyperplane(a,b,C,d);

% check results
res = isemptyobject(hyp1) && ~isemptyobject(hyp2);

% ------------------------------ END OF CODE ------------------------------
