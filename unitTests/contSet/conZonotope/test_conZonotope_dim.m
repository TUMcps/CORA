function res = test_conZonotope_dim
% test_conZonotope_dim - unit test function of dim
%
% Syntax:  
%    res = test_conZonotope_dim
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
% Written:      14-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

% check empty conZonotope
cZ = conZonotope();
if dim(cZ) ~= 0
    res = false;
end


% number of tests
nrOfTests = 1000;

for i=1:nrOfTests
    % random dimension
    n = randi([1,50]);
    % random center
    c = randn(n,1);
    % random generator matrix
    G = 5*randn(n,randi(10));
    nrGens = size(G,2);
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

%------------- END OF CODE --------------