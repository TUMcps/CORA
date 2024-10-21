function res = test_ellipsoid_representsa
% test_ellipsoid_representsa - unit test function of representsa
%
% Syntax:
%    res = test_ellipsoid_representsa
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
% Written:       17-September-2019
% Last update:   ---
% Last revision: 20-July-2023 (MW, rename '...representsa')

% ------------------------------ BEGIN CODE -------------------------------

% 1. comparison to empty set
E = ellipsoid.empty(2);
assert(representsa_(E,'emptySet',eps));
E = ellipsoid([1 0;0 2],[0; 1]);
assert(~representsa_(E,'emptySet',eps));


% 2. comparison to origin
% empty case
E = ellipsoid.empty(2);
assert(~representsa(E,'origin'));

% only origin
E = ellipsoid(zeros(3),zeros(3,1));
assert(representsa(E,'origin'));

% shifted center
E = ellipsoid(zeros(3),0.01*ones(3,1));
assert(~representsa(E,'origin'));

% shifted center, contains origin within tolerance
E = ellipsoid(0.01*eye(3),0.01*ones(3,1));
tol = 0.15;
assert(representsa(E,'origin',tol));


% 3. comparison to capsule
% q = [3;1;-2;1];
% E = ellipsoid(diag(ones(4,1)),q);
% [res(end+1,1),C] = representsa(E,'capsule');
% res(end+1,1) = isequal(C,capsule(q,zeros(4,1),1));


% 4. comparison to point
E = ellipsoid(zeros(4),[3;2;-1;4]);
assert(representsa(E,'point'));
% degenerate
E = ellipsoid([1 0 0; 0 0 0; 0 0 0],[1;2;-1]);
assert(~representsa(E,'point'));
% full-dimensional
E = ellipsoid(eye(4),[3;2;-1;4]);
assert(~representsa(E,'point'));


% 5. comparison to zonotope
E = ellipsoid(zeros(2),[2;1]);
[isZonotope,Z] = representsa(E,'zonotope');
assert(isZonotope)
assert(isequal(Z,zonotope([2;1])));


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
