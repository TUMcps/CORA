function res = testLong_ellipsoid_contains
% testLong_ellipsoid_contains - unit test function of contains
%
% Syntax:
%    res = testLong_ellipsoid_contains
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

% Authors:       Victor Gassmann
% Written:       15-October-2019
% Last update:   07-August-2020
%                19-March-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check ellipsoid 
res = testLong_component_ellipsoid_inEllipsoid;

% check zonotope
res = res && testLong_component_ellipsoid_inZonotope;

% check point containment
runs = 10;
bools = [false,true];
for i=1:5:30
    for j=1:runs
        for k=1:2
            %%% generate all variables necessary to replicate results
            N = 10*i;
            E = ellipsoid.generateRandom('Dimension',i,'IsDegenerate',bools(k));
            samples = randPoint(E,N);
            %%%
            
            if ~all(contains(E,samples))
                res = false; return
            end
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
