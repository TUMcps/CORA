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
rng(1)

% test multiple generateRandom calls;
for i=1:10
    S = contSet.generateRandom();
    assert(isa(S, 'contSet'));
end

% test dimension
S = contSet.generateRandom('Dimension', 3);
assert(dim(S) == 3);

% test given classes
S = contSet.generateRandom({'interval'});
assert(isa(S, 'interval'));

S = contSet.generateRandom({'polyZonotope'});
assert(isa(S, 'polyZonotope'));

S = contSet.generateRandom({'interval','zonotope'});
assert(isa(S, 'interval') || isa(S, 'zonotope'));

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
