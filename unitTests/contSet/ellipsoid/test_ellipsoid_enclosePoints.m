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
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      26-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
res = true;

n = 5;
N = 100;
V = randn(n,N);
[U,S,W] = svd(V);
s = diag(S);
TOL = ellipsoid().TOL;
s(s<=TOL) = 10*TOL;
S = [diag(s),zeros(n,N-n)];
V = U*S*W';
    

%
E = ellipsoid.enclosePoints(V);

% degenerate
s = diag(S);
s(randperm(n,randi([1,n]))) = 0;
Vd = U*[diag(s),zeros(n,N-n)]*W';
Ed = ellipsoid.enclosePoints(Vd);

if ~in(E,V) || ~in(Ed,Vd)
    res = false;
end


if res
    disp([mfilename,' successful']);
else
    disp([mfilename,' failed']);
end
%------------- END OF CODE --------------
