function res = testLong_capsule_isequal
% testLong_capsule_isequal - unit test function of isequal
%
% Syntax:
%    res = testLong_capsule_isequal
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
% Last update:   12-March-2021 (MW, more thorough testing)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

tol = 1e-9;
nrOfTests = 1000;

for i=1:nrOfTests

    % random dimension
    n = randi([2,50]);

    % random properties
    c1 = randn(n,1);
    c2 = randn(n,1);
    g1 = randn(n,1);
    g2 = randn(n,1);
    r1 = rand(1);
    r2 = rand(1);
    C = capsule(c1,g1,r1);

    % test combinations of properties
    % ... different center
    C_ = capsule(c2,g1,r1);
    assertLoop(~isequal(C,C_,tol),i)

    % ... different generator
    C_ = capsule(c1,g2,r1);
    assertLoop(~isequal(C,C_,tol),i)

    % ... different radius
    C_ = capsule(c1,g1,r2);
    assertLoop(~isequal(C,C_,tol),i)

    % ... empty capsule
    C_empty = capsule.empty(n);
    assertLoop(~isequal(C,C_empty),i)

    % ... capsule of reduced dimension
    if n > 1
        C_red = capsule(c1(1:n-1),g1(1:n-1),r1);
        assertLoop(~isequal(C,C_red),i)
    end

end

% ------------------------------ END OF CODE ------------------------------
