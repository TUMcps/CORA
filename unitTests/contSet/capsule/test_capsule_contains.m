function res = test_capsule_contains
% test_capsule_contains - unit test function of contains
%
% Syntax:
%    res = test_capsule_contains
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

% Authors:       Mark Wetzlinger, Adrian Kulmburg
% Written:       27-July-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 1. Same center, generator, different radius
c = [2; -1; 1];
g = [4; -3; 2];
r_small = 0.5; r_big = 1;

% instantiate capsules
C_small = capsule(c,g,r_small);
C_big = capsule(c,g,r_big);

% compute containment
res(end+1,1) = ~contains(C_small,C_big);
res(end+1,1) = contains(C_big,C_small);


% 2. centers too far away from one another
c_plus = 10*c; c_minus = -c;
% norm of generator = 1, radius <= 1
g = [4; -3; 2]; g = g ./ vecnorm(g,2);
r = 0.5;
C_plus = capsule(c_plus,g,r);
C_minus = capsule(c_minus,g,r);

% compute containment
res(end+1,1) = ~contains(C_plus,C_minus);
res(end+1,1) = ~contains(C_minus,C_plus);


% 3. capsules overlapping, different generators
g1 = [4; -3; 2];
g2 = g1 + [-1; 0.2; 0.5];
C1 = capsule(c,g1,r);
C2 = capsule(c,g2,r);
% compute containment
res(end+1,1) = ~contains(C1,C2);
res(end+1,1) = ~contains(C2,C1);


% 4. point containment
% instantiate capsule
C = capsule([1; 1; 0], [0.5; -1; 1], 0.5);
% create a point inside (center)
p_inside = center(C);
% ... and outside
p_outside = 10*ones(dim(C),1);
% check if correct results for containment
res(end+1,1) = contains(C, p_inside);
res(end+1,1) = ~contains(C, p_outside);

% array of points (all outside)
num = 10;
p_array = 10*(ones(dim(C),num)+rand(dim(C),num));
res(end+1,1) = ~any(contains(C, p_array));


% combine tests
res = all(res);


% dimension mismatch
C1 = capsule(1,1,1);
C2 = capsule(rand(2,1),rand(2,1),1);
try
    contains(C1,C2);
catch ME
    if ~strcmp(ME.identifier,'CORA:dimensionMismatch')
        throw(CORAerror('CORA:testFailed'));
    end
end

% ------------------------------ END OF CODE ------------------------------
