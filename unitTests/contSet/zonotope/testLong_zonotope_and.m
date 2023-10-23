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
% Last update:   ---
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
        Z2 = zonotope(ones(dims(d),1),-1+2*rand(dims(d),nrOfGens));
        Z3 = zonotope(100*ones(dims(d),1),-1+2*rand(dims(d),nrOfGens));

        % compute intersection
        Znonempty = Z1 & Z2; % non-empty
        Zempty = Z1 & Z3; % empty
        if representsa(Znonempty,'emptySet') ...
                || ~representsa(Zempty,'emptySet')
            res = false;
            return
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
