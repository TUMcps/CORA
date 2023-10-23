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
% See also: -

% Authors:       Niklas Kochdumper
% Written:       13-January-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% generate random 2D polynomial zonotope
pZ = polyZonotope.generateRandom('Dimension',2,'NrGenerators',10,'NrFactors',4);

% convert to polygon
pgon = polygon(pZ);

% compute interval enclosure
I = interval(pgon);

% generate random points outside the polynomial zonotope and check if "in"
% gives the correct result
for i = 1:100
   
    p = randPoint(I);
    
    if ~contains(pgon,p) && contains(pZ,p,'approx')
        res = false; break
    end
end

% ------------------------------ END OF CODE ------------------------------
