function res = testLong_capsule_isIntersecting
% testLong_capsule_isIntersecting - unit test function of isIntersecting
%    note: only capsule-to-capsule tested
%
% Syntax:
%    res = testLong_capsule_isIntersecting
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
% Written:       11-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

nrTests = 100;

% 1. Same center, generator, different radius

for i=1:nrTests

    % random dimension
    n = randi([2,10]);

    c = randn(n,1);
    g = randomPointOnSphere(n);
    r_small = rand(1);
    r_big = 2*r_small;

    % instantiate capsules
    C_small = capsule(c,g,r_small);
    C_big = capsule(c,g,r_big);

    % check for correctness
    assert(isIntersecting(C_small,C_big));
    assert(isIntersecting(C_big,C_small));

end


% 2. centers too far away from one another

for i=1:nrTests

    % random dimension
    n = randi([2,10]);

    % center far away
    c_plus = 10*rand(n,1)+10;
    c_minus = -c_plus;
    % norm of generator = 1, radius <= 1
    g = randomPointOnSphere(n);
    r = rand(1);

    % instantiate capsules
    C_plus = capsule(c_plus,g,r);
    C_minus = capsule(c_minus,g,r);

    % check for correctness
    assert(~isIntersecting(C_plus,C_minus));
    assert(~isIntersecting(C_minus,C_plus));

end


% 3. capsules overlapping

for i=1:nrTests

    % random dimension
    n = randi([2,10]);
    
    % same center and radius
    c = randn(n,1);
    r = rand(1);
    % different generators
    g1 = randomPointOnSphere(n);
    g2 = randomPointOnSphere(n);

    % instantiate capsules
    C1 = capsule(c,g1,r);
    C2 = capsule(c,g2,r);

    % check for correctness
    assert(isIntersecting(C1,C2));
    assert(isIntersecting(C2,C1));

end


% 4. capsules touching in exactly one point

for i=1:nrTests
    
    % random dimension
    n = randi([2,10]);

    % same radius
    r = randi(5);
    % random generator
    g1 = randomPointOnSphere(n);
    g2 = -g1;
    % set centers so that intersection is the origin
    c1 = (1+r)*g1;
    c2 = (1+r)*g2;

    % instantiate capsules
    C1 = capsule(c1,g1,r);
    C2 = capsule(c2,g2,r);

    % check for correctness
    assert(isIntersecting(C1,C2));
    assert(isIntersecting(C2,C1));

end

% combine tests
res = true;

% ------------------------------ END OF CODE ------------------------------
