function res = testLongDuration_zonotope_quadMap
% testLongDuration_zonotope_quadMap - unit test function of quadMap
%
% Syntax:  
%    res = testLongDuration_zonotope_quadMap
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Matthias Althoff, Mark Wetzlinger
% Written:      26-July-2016
% Last update:  09-August-2020 (MW, extend by random points)
% Last revision:---

%------------- BEGIN CODE --------------

% create zonotopes
Z1 = zonotope([-4, -3, -2; 1, 2, 3]);
Z2 = zonotope([1, 4, 2; -3, 2, -1]);

nrOfRandPoints = 1000;

% create matrices
Q{1} = [1 -2; -3 4];
Q{2} = [0.5 0; 2 -1];

Zres = quadMap(Z1,Q);

% map random points in zonotope and check if they are inside the result
res(1) = true;
for i=1:nrOfRandPoints
    p = randPoint(Z1);
    for d=1:dim(Z1)
        pQp(d,1) = p' * Q{d} * p;
    end
    if ~in(Zres,pQp)
        res(1) = false;
        break
    end
end

% 2. quadMapMixed: Z1*Q*Z2 ------------------------------------------------

% obtain result
Zres = quadMap(Z1,Z2,Q);

% map random points in zonotope and check if they are inside the result
res(2) = true;
for i=1:nrOfRandPoints
    p1 = randPoint(Z1);
    p2 = randPoint(Z2);
    for d=1:dim(Z1)
        pQp(d,1) = p1' * Q{d} * p2;
    end
    if ~in(Zres,pQp)
        res(2) = false;
        break
    end
end


% gather results
res = all(res);

if res
    disp('testLongDuration_zonotope_quadMap successful');
else
    disp('testLongDuration_zonotope_quadMap failed');
end

%------------- END OF CODE --------------
