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
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Mark Wetzlinger
% Written:      27-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

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

% compute intersection
bigInSmall = isIntersecting(C_small,C_big);
smallInBig = isIntersecting(C_big,C_small);

% check for correctness
if ~bigInSmall || ~smallInBig
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

% compute intersection
minusInPlus = isIntersecting(C_plus,C_minus);
plusInMinus = isIntersecting(C_minus,C_plus);

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

% compute intersection
C2inC1 = isIntersecting(C1,C2);
C1inC2 = isIntersecting(C2,C1);

% check for correctness
if ~C2inC1 || ~C1inC2
    res(3) = false;
end


% 4. capsules touching in exactly one point
res(4) = true;

% same radius
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

% compute intersection
C2inC1 = isIntersecting(C1,C2);
C1inC2 = isIntersecting(C2,C1);

% check for correctness
if ~C2inC1 || ~C1inC2
    res(4) = false;
end
    

% 5. dimension mismatch
res(5) = true;
C1 = capsule(1,1,1);
C2 = capsule(rand(2,1),rand(2,1),1);
try
    isIntersecting(C1,C2);
catch ME
    if ~strcmp(ME.identifier,'CORA:dimensionMismatch')
        res(5) = false;
    end
end


% combine tests
res = all(res);

if res
    disp('test_isIntersecting successful');
else
    disp('test_isIntersecting failed');
end

%------------- END OF CODE --------------