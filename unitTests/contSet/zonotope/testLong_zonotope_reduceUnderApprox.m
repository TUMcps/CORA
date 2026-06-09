function res = testLong_zonotope_reduceUnderApprox
% testLong_zonotope_reduceUnderApprox - unit test of reduceUnderApprox
%
% Syntax:
%    res = testLong_zonotope_reduceUnderApprox
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

% Authors:       Niklas Kochdumper
% Written:       17-January-2025
% Last update:   10-February-2026 (LL, new methods & tol for contains)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% number of tests
nrOfTests = 3;

% order reduction methods
methods = {'sadraddini','yang','raghuraman','kochdumper',...
                    'scale','boxCone','nlp','cluster'};

% loop over all test cases
for i = 1:nrOfTests

    % random dimension
    n = randi([2,5]); % small because of containment checks

    % random original and reduced zonotope order
    order = randi([2,5]);
    orderRed = randi([1,order-1]);

    % random zonotope
    Z = zonotope.generateRandom("Dimension",n,"NrGenerators",order*n);

    comb = combinator(size(Z.G,2),n,'c');
    
    % loop over all order reduction methods
    for j = 1:length(methods)

        % compute reduced order zonotope
        Zred = reduceUnderApprox(Z,methods{j},orderRed);

        % check that the number of generators is correct
        resGen = size(Zred.G,2) == orderRed*n;
        assertLoop(resGen,i,j)

        % downscale reduced zonotope minimally for numerical robustness of
        % containment checks
        Zred = zonotope(center(Zred), 0.9999*generators(Zred));

        % check that the reduced order zonotope is contained inside the
        % original zonotope        
        if size(comb,1) <= 1000
            resCon = contains(polytope(Z),Zred);
        else
            p = randPoint(Zred,100,'extreme');
            resCon = all(contains(Z,p,'exact',1e-4));
        end

        % check for containment in zonotope
        assertLoop(resCon,i,j)
    end
end

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
