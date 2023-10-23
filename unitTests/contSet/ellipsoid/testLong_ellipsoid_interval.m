function res = testLong_ellipsoid_interval
% testLong_ellipsoid_interval - unit test function of interval
%
% Syntax:
%    res = testLong_ellipsoid_interval
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
% Written:       16-October-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
nRuns = 2;

for i=10:5:15
    for j=1:nRuns
        E = ellipsoid.generateRandom('Dimension',i,'IsDegenerate',false);
        Y = randPoint(E,1000);
        %compute interval
        I = interval(E);
        %check if all points are in the interval
        if ~contains(I,Y)
            res = false; return
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
