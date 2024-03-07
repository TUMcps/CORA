function res = testLong_zonotope_and
% testLong_zonotope_and - unit test function of and
%
% Syntax:
%    res = testLong_zonotope_and
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
% Written:       09-September-2020
% Last update:    03-March-2024 (TL, bug fix, Z2 now always intersects with Z1)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% dimensions
dims = 2:2:8;
% number of tests
testsPerDim = 10;

% box has to be the same as conversion to interval
for d=1:length(dims)
    for test=1:testsPerDim
        % create random zonotopes
        nrOfGens = 10;
        Z1 = zonotope(zeros(dims(d),1),-1+2*rand(dims(d),nrOfGens)); 

        % generate second zonotope on boundary (has to intersect)
        x_boundary = Z1.randPoint(1,'boundary');
        Z2 = zonotope(x_boundary,-1+2*rand(dims(d),nrOfGens));

        % generate third zonotope so far away that they cannot intersect
        Z3 = zonotope(100*ones(dims(d),1),-1+2*rand(dims(d),nrOfGens));

        % compute intersection
        Znonempty = Z1 & Z2; % non-empty
        Zempty = Z1 & Z3; % empty

        % test assumptions
        if representsa(Znonempty,'emptySet') ...
                || ~representsa(Zempty,'emptySet')
            res = false;
            return
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
