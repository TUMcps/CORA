function res = test_conHyperplane_generateRandom
% test_conHyperplane_generateRandom - unit test function of generateRandom
%
% Syntax:
%    res = test_conHyperplane_generateRandom
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
% Written:       02-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

% empty case
hyp = conHyperplane.generateRandom();
resvec(end+1) = isa(hyp,'conHyperplane');

% dimension given
hyp = conHyperplane.generateRandom('Dimension',3);
resvec(end+1) = dim(hyp) == 3;

% normal vector given
hyp = conHyperplane.generateRandom('NormalVector',[1; 2]);
resvec(end+1) = all(hyp.a == [1 2]);

% offset given
hyp = conHyperplane.generateRandom('Offset',3);
resvec(end+1) = all(hyp.b == 3);

% unify results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
