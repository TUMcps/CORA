function res = testLong_capsule_contains
% testLong_capsule_contains - unit test function of contains
%
% Syntax:
%    res = testLong_capsule_contains
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

res = true;
nrTests = 1000;

% 1. Same center, generator, different radius

for i=1:nrTests

    % random dimension
    n = randi([2,30]);

    c = randn(n,1);
    g = randomPointOnSphere(n);
    r_small = rand(1);
    r_big = 2*r_small;

    % instantiate capsules
    C_small = capsule(c,g,r_small);
    C_big = capsule(c,g,r_big);

    % check for correctness
    assert(~contains(C_small,C_big));
    assert(contains(C_big,C_small));

end


% 2. centers too far away from one another

for i=1:nrTests

    % random dimension
    n = randi([2,30]);

    % center far away
    c_plus = 10*rand(n,1);
    c_minus = -c_plus;
    % norm of generator = 1, radius <= 1
    g = randomPointOnSphere(n);
    r = rand(1);

    % instantiate capsules
    C_plus = capsule(c_plus,g,r);
    C_minus = capsule(c_minus,g,r);

    % check correctness
    assert(~contains(C_plus,C_minus));
    assert(~contains(C_minus,C_plus));

end


% 3. capsules overlapping

for i=1:nrTests

    % random dimension
    n = randi([2,30]);

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
    assert(~contains(C1,C2));
    assert(~contains(C2,C1));
end

% ------------------------------ END OF CODE ------------------------------
