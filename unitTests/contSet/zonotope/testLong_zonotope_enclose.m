function res = testLong_zonotope_enclose
% testLong_zonotope_enclose - unit test function of enclose
%
% Syntax:
%    res = testLong_zonotope_enclose
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

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       26-July-2016
% Last update:   09-August-2020 (MW, enhance randomness)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% number of tests
nrTests = 100;
ptsPerLine = 10;

for i=1:nrTests

    % random dimension
    n = randi([2,4]);

    % create two random zonotopes
    nrOfGens = randi([5,15],1,1);
    c1 = -20*ones(n,1);
    G1 = -1+2*rand(n,nrOfGens);
    Z1 = zonotope(c1,G1);
    c2 = 20*ones(n,1);
    G2 = -1+2*rand(n,nrOfGens);
    Z2 = zonotope(c2,G2);
    
    % compute enclosure
    Zenc = enclose(Z1,Z2);

    % random points in Z1 or Z2
    p1 = randPoint(Z1);
    p2 = randPoint(Z2);
    % connect points by line, all have to be in enclosure
    pts = p1 + (p2-p1) .* linspace(0,1,ptsPerLine);

    % random points have to be also in Zenc
    if ~all(contains(Zenc,pts))
        res = false; return
    end
end

% ------------------------------ END OF CODE ------------------------------
