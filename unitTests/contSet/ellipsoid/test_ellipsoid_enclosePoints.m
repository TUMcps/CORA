function res = test_ellipsoid_enclosePoints
% test_ellipsoid_enclosePoints - unit test function of enclosePoints
%
% Syntax:
%    res = test_ellipsoid_enclosePoints
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
% Written:       26-July-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% dimension
n = 5;
% number of points
N = 100;
V = randn(n,N);

% svd
[U,S,W] = svd(V);
s = diag(S);

TOL = 1e-6; % default value from constructor
s(s<=TOL) = 10*TOL;
S = [diag(s),zeros(n,N-n)];
V = U*S*W';
    

% enclose points
E = ellipsoid.enclosePoints(V);

% degenerate (n-1)
s = diag(S);
s(randperm(n,randi([1,n]))) = 0;
Vd = U*[diag(s),zeros(n,N-n)]*W';
Ed = ellipsoid.enclosePoints(Vd);

if ~all(contains(E,V)) || ~all(contains(Ed,Vd))
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
