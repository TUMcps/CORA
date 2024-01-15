function res = test_contSet_generateRandom
% test_contSet_generateRandom - unit test function of
%    contSet.generateRandom
%
% Syntax:
%    res = test_contSet_generateRandom
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
% Written:       05-April-2023
% Last update:   ---
% Last revision: 09-January-2024

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true(0);
rng(1)

% test multiple generateRandom calls;
for i=1:10
    S = contSet.generateRandom();
    res(end+1) = isa(S, 'contSet');
end

% test dimension
S = contSet.generateRandom('Dimension', 3);
res(end+1) = dim(S) == 3;

% test given classes
S = contSet.generateRandom({'interval'});
res(end+1) = isa(S, 'interval');

S = contSet.generateRandom({'polyZonotope'});
res(end+1) = isa(S, 'polyZonotope');

S = contSet.generateRandom({'interval','zonotope'});
res(end+1) = isa(S, 'interval') || isa(S, 'zonotope');

% test if all were successfull
res = all(res);

% ------------------------------ END OF CODE ------------------------------
