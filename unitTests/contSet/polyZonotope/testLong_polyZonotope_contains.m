function res = testLong_polyZonotope_contains
% testLong_polyZonotope_contains - unit test function for containment
%    checks of polynomial zonotopes
%
% Syntax:
%    res = testLong_polyZonotope_contains
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
% Written:       13-January-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% generate random 2D polynomial zonotope
pZ = polyZonotope.generateRandom('Dimension',2,'NrGenerators',10,'NrFactors',4);

% enclose the polynomial zonotope by a polygon
pgon = polygon(pZ);

% enclose the polygon by an interval
I = interval(pgon);

% sample points within the interval
p = randPoint(I,100);

% for each randomly sampled point, if it is not contained in the polygon, 
% it must not be contained in the polynomial zonotope
for i = 1:size(p,2)
    if ~contains(pgon,p(:,i))
        assert(~contains(pZ,p(:,i),'approx'))
    end
end

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
