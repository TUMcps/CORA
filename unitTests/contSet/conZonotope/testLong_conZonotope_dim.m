function res = testLong_conZonotope_dim
% testLong_conZonotope_dim - unit test function of dim
%
% Syntax:
%    res = testLong_conZonotope_dim
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
% Written:       03-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% number of tests
nrOfTests = 1000;

for i=1:nrOfTests
    
    % random dimension
    n = randi(50);

    % random number of generators
    nrGens = randi(2*n);

    % random center and generator matrix
    c = randn(n,1);
    G = 5*randn(n,nrGens);
    
    % random constraints
    A = diag(randn(nrGens,1));
    b = randn(nrGens,1);
    
    % instantiate conZonotope
    cZ = conZonotope(c,G,A,b);
    
    % get dimension
    Zdim = dim(cZ);
    
    % assert correctness
    if Zdim ~= n
        res = false; break
    end
end

% ------------------------------ END OF CODE ------------------------------
