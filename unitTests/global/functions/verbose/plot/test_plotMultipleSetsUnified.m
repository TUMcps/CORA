function res = test_plotMultipleSetsUnified
% test_plotMultipleSetsUnified - unit test function for plotMultipleSetsUnified
%
% Syntax:
%    res = test_plotMultipleSetsUnified
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

% Authors:       Tobias Ladner
% Written:       11-July-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% create sets
N = 20;
sets = arrayfun(@(i) zonotope.generateRandom('Dimension',2),1:N,'UniformOutput',false);
% make sure they are all connected
for i=1:(N-1)
    % sample two points and make sure the subsequent set matches the
    % previous one
    S1 = sets{i}; p1 = randPoint(S1);
    S2 = sets{i+1}; p2 = randPoint(S2);
    % shift
    S2 = S2 + p1-p2;
    sets{i+1} = S2;
end

try
    figure;
    ax = gca(); hold off;

    % plot unified
    numTotalSets = 1;
    plotMultipleSetsUnified(sets,1:2,{'UnifyTotalSets',numTotalSets});
    assert(numel(allchild(ax)) == numTotalSets);

    % plot individually
    numTotalSets = 2;
    plotMultipleSetsUnified(sets,1:2,{'UnifyTotalSets',numTotalSets});
    assert(numel(allchild(ax)) == numTotalSets);

    % plot individually
    numTotalSets = 10;
    plotMultipleSetsUnified(sets,1:2,{'UnifyTotalSets',numTotalSets});
    assert(numel(allchild(ax)) == numTotalSets);

    % plot individually
    numTotalSets = N;
    plotMultipleSetsUnified(sets,1:2,{'UnifyTotalSets',numTotalSets});
    assert(numel(allchild(ax)) == numTotalSets);

    % plot automatically
    numTotalSets = [];
    plotMultipleSetsUnified(sets,1:2,{'UnifyTotalSets',numTotalSets});
    assert(numel(allchild(ax)) < N);

    close;

catch ME
    close
    rethrow(ME)
end

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
