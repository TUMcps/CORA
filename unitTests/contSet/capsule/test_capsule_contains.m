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

% 1. Same center, generator, different radius
res(1) = true;

% define capsule
c = [2; -1; 1];
g = [4; -3; 2];
r_small = 0.5;
r_big = 1;

% instantiate capsules
C_small = capsule(c,g,r_small);
C_big = capsule(c,g,r_big);

% compute containment
bigInSmall = contains(C_small,C_big);
smallInBig = contains(C_big,C_small);

% check for correctness
if bigInSmall || ~smallInBig
    res(1) = false;
end


% 2. centers too far away from one another
res(2) = true;

% center far away
c_plus = 10*c;
c_minus = -c;
% norm of generator = 1, radius <= 1
g = g ./ vecnorm(g,2);
r = 0.5;

% instantiate capsules
C_plus = capsule(c_plus,g,r);
C_minus = capsule(c_minus,g,r);

% compute containment
minusInPlus = contains(C_plus,C_minus);
plusInMinus = contains(C_minus,C_plus);

% check for correctness
if minusInPlus || plusInMinus
    res(2) = false;
end


% 3. capsules overlapping
res(3) = true;

% different generators
g1 = g;
g2 = g + [-1; 0.2; 0.5];

% instantiate capsules
C1 = capsule(c,g1,r);
C2 = capsule(c,g2,r);

% compute containment
C2inC1 = contains(C1,C2);
C1inC2 = contains(C2,C1);

% check for correctness
if C2inC1 || C1inC2
    res(3) = false;
end

% 4. dimension mismatch
res(4) = true;
C1 = capsule(1,1,1);
C2 = capsule(rand(2,1),rand(2,1),1);
try
    contains(C1,C2);
catch ME
    if ~strcmp(ME.identifier,'CORA:dimensionMismatch')
        res(4) = false;
    end
end

% 5. point containment
% instantiate capsule
n = 3;
C = capsule([1; 1; 0], [0.5; -1; 1], 0.5);

% Single points
% create a point inside (center)
p_inside = center(C);
% ... and outside
p_outside = 10*ones(n,1);

% check if correct results for containment
res_inside  = contains(C, p_inside);
res_outside = contains(C, p_outside);

% compare results
res_single = res_inside && ~res_outside;

% Array of points (all outside)
num = 10;
p_array = 10*(ones(n,num)+rand(n,num));
res_array = contains(C, p_array);

res(5) = res_single && ~any(res_array);

% combine tests
res = all(res);

% ------------------------------ END OF CODE ------------------------------
