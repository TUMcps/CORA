function res = test_capsule_representsa
% test_capsule_representsa - unit test function of representsa
%
% Syntax:
%    res = test_capsule_representsa
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
C = capsule.empty(2);
res(end+1,1) = representsa(C,'emptySet');
C = capsule([1;1],[0;1],0.5);
res(end+1,1) = ~representsa(C,'emptySet');
C = capsule([1;1],[0;1],0);
res(end+1,1) = ~representsa(C,'emptySet');


% 2. comparison to interval

% empty capsule
C = capsule.empty(2);
res(end+1,1) = representsa(C,'interval');

% full-dimensional capsule
c = [2; 0; -1];
g = [0.2; -0.7; 0.4];
r = 1;
C = capsule(c,g,r);
res(end+1,1) = ~representsa(C,'interval');

% one-dimensional capsule
C = capsule(2,1,0);
[res(end+1,1),I] = representsa(C,'interval');
res(end+1,1) = isequal(I,interval(1,3));

% two-dimensional capsule with axis-aligned generator and no radius
C = capsule([1;-1],[1;0],0);
[res(end+1,1),I] = representsa(C,'interval');
res(end+1,1) = isequal(I,interval([0;-1],[2;-1]));

% two-dimensional capsule with all-zero generator and (no) radius
C = capsule([0;-1],[0;0],1);
res(end+1,1) = ~representsa(C,'interval');
C = capsule([0;-1],[0;0],0);
[res(end+1,1),I] = representsa(C,'interval');
res(end+1,1) = isequal(I,interval([0;-1]));


% 3. comparison to origin

% empty case
C = capsule.empty(2);
res(end+1,1) = ~representsa(C,'origin');

% true cases
C = capsule(zeros(3,1));
res(end+1) = representsa(C,'origin');
C = capsule(zeros(3,1),zeros(3,1),0);
res(end+1) = representsa(C,'origin');

% shifted center
C = capsule(ones(3,1),zeros(3,1),0);
res(end+1) = ~representsa(C,'origin');

% including generator, no radius
C = capsule(zeros(2,1),ones(2,1),0);
res(end+1) = ~representsa(C,'origin');

% no generator, but radius
C = capsule(zeros(4,1),zeros(4,1),1);
res(end+1) = ~representsa(C,'origin');

% not zero, but within tolerance
C = capsule(zeros(3,1),zeros(3,1),1);
tol = 2;
res(end+1) = representsa(C,'origin',tol);

% does not contain origin, but within tolerance
c = 0.5*[sqrt(2); sqrt(2)];
g = 0.5*[sqrt(2); -sqrt(2)];
r = 0.5;
tol = 2;
C = capsule(c,g,r);
res(end+1) = representsa(C,'origin',tol);


% 4. comparison to ellipsoid
% only center
C = capsule([3;1;2]);
[res(end+1,1),E] = representsa(C,'ellipsoid');
res(end+1,1) = isequal(E,ellipsoid(zeros(3),[3;1;2]));

% ball
C = capsule([4;1;-2;3],zeros(4,1),2);
[res(end+1,1),E] = representsa(C,'ellipsoid');
res(end+1,1) = isequal(E,ellipsoid(4*eye(4),[4;1;-2;3]));

% no ball
C = capsule([3;1],[-1;2],1);
res(end+1,1) = ~representsa(C,'ellipsoid');


% 5. comparison to zonotope
% only center
C = capsule([3;1;2]);
[res(end+1,1),Z] = representsa(C,'zonotope');
res(end+1,1) = isequal(Z,zonotope([3;1;2]));

% center and one generator
C = capsule([3;1;2],[-1;1;0]);
[res(end+1,1),Z] = representsa(C,'zonotope');
res(end+1,1) = isequal(Z,zonotope([3;1;2],[-1;1;0]));

% regular capsule
C = capsule([-1;2],[1;-1],1);
res(end+1,1) = ~representsa(C,'zonotope');


% 6. comparison to point
C = capsule([3;1;2]);
res(end+1,1) = representsa(C,'point');
C = capsule([-1;2],[1;-1],1);
res(end+1,1) = ~representsa(C,'point');


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
