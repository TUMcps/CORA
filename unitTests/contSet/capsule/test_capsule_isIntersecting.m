function res = test_capsule_isIntersecting
% test_capsule_isIntersecting - unit test function of isIntersecting
%    note: only capsule-to-capsule tested
%
% Syntax:
%    res = test_capsule_isIntersecting
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
% Written:       27-July-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% empty capsule
C_empty = capsule.empty(2);
c = [2;1]; g = [0;-1]; r = 2;
C = capsule(c,g,r);
res(end+1,1) = ~isIntersecting(C_empty,C_empty);
res(end+1,1) = ~isIntersecting(C_empty,C);


% same center, generator, different radius
c = [2; -1; 1]; g = [4; -3; 2];
r_small = 0.5; r_big = 1;
% instantiate capsules
C_small = capsule(c,g,r_small);
C_big = capsule(c,g,r_big);
% check intersection
res(end+1,1) = isIntersecting(C_small,C_big);
res(end+1,1) = isIntersecting(C_big,C_small);


% centers too far away from one another
c_plus = 10*c;
c_minus = -c;
% norm of generator = 1, radius <= 1
g = [4; -3; 2]; g = g ./ vecnorm(g,2);
r = 0.5;
% instantiate capsules
C_plus = capsule(c_plus,g,r);
C_minus = capsule(c_minus,g,r);
% check intersection
res(end+1,1) = ~isIntersecting(C_plus,C_minus);
res(end+1,1) = ~isIntersecting(C_minus,C_plus);


% overlapping capsules
c = [2; -1; 1]; r = 0.5;
g1 = g;
g2 = g + [-1; 0.2; 0.5];
% instantiate capsules
C1 = capsule(c,g1,r);
C2 = capsule(c,g2,r);
% check intersection
res(end+1,1) = isIntersecting(C1,C2);
res(end+1,1) = isIntersecting(C2,C1);


% capsules touching in exactly one point, same radius
r = 3;
% random generator (unit length for simplicity)
g1 = [3; 4] / 5;
g2 = -g1;
% set centers so that intersection is the origin
c1 = (1+r)*g1;
c2 = (1+r)*g2;
% instantiate capsules
C1 = capsule(c1,g1,r);
C2 = capsule(c2,g2,r);
% check intersection
res(end+1,1) = isIntersecting(C1,C2);
res(end+1,1) = isIntersecting(C2,C1);


% combine tests
res = all(res);


% dimension mismatch
C1 = capsule(1,1,1);
C2 = capsule(rand(2,1),rand(2,1),1);
try
    isIntersecting(C1,C2);
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
