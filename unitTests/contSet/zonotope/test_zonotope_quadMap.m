function res = test_zonotope_quadMap
% test_zonotope_quadMap - unit test function of
%  quadMap
%
% Syntax:  
%    res = test_zonotope_quadMap
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

% create matrices
Q{1} = [1 -2; -3 4];
Q{2} = [0.5 0; 2 -1];

% random testing
nrOfRandPoints = 5000;

% 1. quadMapSingle: Z1*Q*Z1 -----------------------------------------------

% obtain result
Zres = quadMap(Z1,Q);

% obtain zonotope matrix
Zmat = Zres.Z;

% true result
true_mat = [102.5, 27.5, 35, 95, 110, 125; ...
            -16.25, -5.75, -9.5, -14, -26, -32];

% compare solutions
res_comp(1) = all(all(Zmat == true_mat));

% map random points in zonotope and check if they are inside the result
res_randpoint(1) = true;
for i=1:nrOfRandPoints
    p = randPoint(Z1);
    for d=1:dim(Z1)
        pQp(d,1) = p' * Q{d} * p;
    end
    if ~containsPoint(Zres,pQp)
        res_randpoint(1) = false;
        break
    end
end

% 2. quadMapMixed: Z1*Q*Z2 ------------------------------------------------

% obtain result
Zres = quadMap(Z1,Z2,Q);

% obtain zonotope matrix
Zmat = Zres.Z;

% true result
true_mat = [-43, -51, -59, -4, -8, -12, -26, -32, -38;
            3, 8.5, 14, -2, 6, 14, 1, 7, 13];

% compare solutions
res_comp(2) = all(all(Zmat == true_mat));

% map random points in zonotope and check if they are inside the result
res_randpoint(2) = true;
for i=1:nrOfRandPoints
    p1 = randPoint(Z1);
    p2 = randPoint(Z2);
    for d=1:dim(Z1)
        pQp(d,1) = p1' * Q{d} * p2;
    end
    if ~containsPoint(Zres,pQp)
        res_randpoint(2) = false;
        break
    end
end


% gather results
res = all(res_comp) && all(res_randpoint);

if res
    disp('test_zonotope_quadMap successful');
else
    disp('test_zonotope_quadMap failed');
end

%------------- END OF CODE --------------
