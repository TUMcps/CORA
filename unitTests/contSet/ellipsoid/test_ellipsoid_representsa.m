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

res = true(0);

% 1. comparison to empty set
E = ellipsoid.empty(2);
res(end+1,1) = representsa_(E,'emptySet',eps);
E = ellipsoid([1 0;0 2],[0; 1]);
res(end+1,1) = ~representsa_(E,'emptySet',eps);


% 2. comparison to origin
% empty case
E = ellipsoid.empty(2);
res(end+1,1) = ~representsa(E,'origin');

% only origin
E = ellipsoid(zeros(3),zeros(3,1));
res(end+1,1) = representsa(E,'origin');

% shifted center
E = ellipsoid(zeros(3),0.01*ones(3,1));
res(end+1,1) = ~representsa(E,'origin');

% shifted center, contains origin within tolerance
E = ellipsoid(0.01*eye(3),0.01*ones(3,1));
tol = 0.15;
res(end+1,1) = representsa(E,'origin',tol);


% 3. comparison to capsule
% q = [3;1;-2;1];
% E = ellipsoid(diag(ones(4,1)),q);
% [res(end+1,1),C] = representsa(E,'capsule');
% res(end+1,1) = isequal(C,capsule(q,zeros(4,1),1));


% 4. comparison to point
E = ellipsoid(zeros(4),[3;2;-1;4]);
res(end+1,1) = representsa(E,'point');


% 5. comparison to zonotope
E = ellipsoid(zeros(2),[2;1]);
[res(end+1,1),Z] = representsa(E,'zonotope');
res(end+1,1) = isequal(Z,zonotope([2;1]));


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
