function res = testLong_zonotope_quadMap
% testLong_zonotope_quadMap - unit test function of quadMap
%
% Syntax:
%    res = testLong_zonotope_quadMap
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
% Last update:   09-August-2020 (MW, extend by random points)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% create zonotopes
Z1 = zonotope([-4, -3, -2; 1, 2, 3]);
Z2 = zonotope([1, 4, 2; -3, 2, -1]);

nrOfRandPoints = 1000;

% create matrices
Q{1} = [1 -2; -3 4];
Q{2} = [0.5 0; 2 -1];

Zres = quadMap(Z1,Q);

% map random points in zonotope and check if they are inside the result
for i=1:nrOfRandPoints
    p = randPoint(Z1);
    pQp = quadMapPoint(p,p,Q);
    if ~contains(Zres,pQp)
        res = false; return
    end
end

% 2. quadMapMixed: Z1*Q*Z2 ------------------------------------------------

% obtain result
Zres = quadMap(Z1,Z2,Q);

% map random points in zonotope and check if they are inside the result
for i=1:nrOfRandPoints
    p1 = randPoint(Z1);
    p2 = randPoint(Z2);
    pQp = quadMapPoint(p1,p2,Q);
    if ~contains(Zres,pQp)
        res = false; return
    end
end

% ------------------------------ END OF CODE ------------------------------
