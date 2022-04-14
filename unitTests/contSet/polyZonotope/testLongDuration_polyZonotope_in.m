function res = testLongDuration_polyZonotope_in
% testLongDuration_polyZonotope_in - unit test function for containment
%    checks of polynomial zonotopes
%
% Syntax:  
%    res = testLongDuration_polyZonotope_in
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Niklas Kochdumper
% Written:      13-January-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = false;


% Random Test -------------------------------------------------------------

% generate random 2D polynomial zonotope
pZ = polyZonotope.generateRandom(2,10,4);

% convert to polygon
pgon = polygon(pZ);

% compute interval enclosure
int = interval(pgon);

% generate random points outside the polynomial zonotope and check if "in"
% gives the correct result
for i = 1:100
   
    p = randPoint(int);
    
    if ~in(pgon,p)
        if in(pZ,p,'approx')
            error('testLongDuration_polyZonotope_in: random test failed!');
        end
    end
end

res = true;

%------------- END OF CODE --------------