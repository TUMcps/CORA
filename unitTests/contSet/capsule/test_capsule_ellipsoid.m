function res = test_capsule_ellipsoid
% test_capsule_ellipsoid - unit test function of conversion to ellipsoids
%
% Syntax:
%    res = test_capsule_ellipsoid
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

% Authors:       Mark Wetzlinger
% Written:       25-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty set case
n = 2;
C = capsule.empty(n);
% convert to ellipsoid
E = ellipsoid(C);
res = representsa_(E,'emptySet',eps) && dim(E) == 2;

% ball case
C = capsule([1;-1],[0;0],3);
E = ellipsoid(C);
% check
res(end+1,1) = isequal(E,ellipsoid(9*eye(2),[1;-1]));

% 2D capsule
C = capsule([1;-1],[6;2],2);
% conversion to outer-/inner-approximation
Eo = ellipsoid(C,'outer');
Ei = ellipsoid(C,'inner');
% check containment
res(end+1,1) = contains(Eo,C);
% res(end+1,1) = contains(C,Ei); % method not precise enough...

% sample points and check for containment
% pC = randPoint(C,100,'standard');
pEi = randPoint(Ei,100,'standard');
% res(end+1,1) = all(contains(Eo,pC));
res(end+1,1) = all(contains(C,pEi));

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
