function res = testLong_capsule_capsule
% testLong_capsule_capsule - unit test function of capsule (constructor)
%
% Syntax:
%    res = testLong_capsule_capsule
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
% See also: none

% Authors:       Mark Wetzlinger
% Written:       19-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
nrOfTests = 50;
for i=1:nrOfTests

    % random dimension
    n = randi([2,25]);
    
    % random center, generator, and radius
    c = randn(n,1);
    g = randn(n,1);
    r = rand(1);
    
    % admissible initializations
    % only center
    C = capsule(c);
    assertLoop(compareMatrices(C.c,c),i)

    % center and generator
    C = capsule(c,g);
    assertLoop(compareMatrices(C.c,c),i)
    assertLoop(compareMatrices(C.g,g),i)
    
    % center, generator, and radius
    C = capsule(c,g,r);
    assertLoop(compareMatrices(C.c,c),i)
    assertLoop(compareMatrices(C.g,g),i)
    assertLoop(withinTol(C.r,r),i)
    
    % prevent 1D cases (where g can be interpreted as r)
    n = n + 1;
    c = randn(n,1);
    g = randn(n,1);
    
    % wrong initializations
    cn_1 = randn(n+1,1);
    gn_1 = randn(n+1,1);
    rneg = -rand(1);
    rvec = rand(n,1);
    
    % mismatch between center and generator
    assertThrowsAs(@capsule,'CORA:wrongInputInConstructor',c,gn_1);
    assertThrowsAs(@capsule,'CORA:wrongInputInConstructor',cn_1,g);
        
    % negative radius
    assertThrowsAs(@capsule,'CORA:wrongValue',c,g,rneg);
    
    % radius as a vector
    assertThrowsAs(@capsule,'CORA:wrongValue',c,g,rvec);
    
    % too many input arguments
    assertThrowsAs(@capsule,'CORA:numInputArgsConstructor',c,g,r,r);
    
end

% ------------------------------ END OF CODE ------------------------------
