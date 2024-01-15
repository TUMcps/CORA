function res = test_capsule_enclosePoints
% test_capsule_enclosePoints - unit test function of enclosePoints
%
% Syntax:
%    res = test_capsule_enclosePoints
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

% Authors:       Mark Wetzlinger
% Written:       23-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty set case
n = 2;
p = zeros(n,0);
C = capsule.enclosePoints(p);
res = representsa_(C,'emptySet',eps) && dim(C) == n;

% points
p = [-1  1 2 3 2  1 -2 -4 3;...
      2 -1 3 4 1 -1  2 -3 2];
% compute enclosing capsule
C = capsule.enclosePoints(p);

% check if all points are contained in capsule
res(end+1,1) = all(contains(C,p));

% only one point
p = [-1;1;3];
C = capsule.enclosePoints(p);
% center has to be point, generator and radius 0
res(end+1,1) = all(p == center(C)) && C.r == 0 && ~any(C.g);

% visualization
% figure; hold on;
% plot(C);
% plot(p(1,:),p(2,:),'.k');

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
