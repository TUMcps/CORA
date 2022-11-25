function res = testLongDuration_conZonotope_isempty
% testLongDuration_conZonotope_isempty - unit test function of isempty
%
% Syntax:  
%    res = testLongDuration_conZonotope_isempty
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

% check empty conZonotope object
conZono = conZonotope();
if ~isempty(conZono)
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
    G = randn(n,n+randi(10));
    nrGens = size(G,2);
    
    % instantiate conZonotope without constraints
    conZono = conZonotope(c,G);
    
    % assert correctness
    if isempty(conZono)
        res = false; break;
    end
    
    % random constraints so that conZonotope represents just a point
    % as A being diagional forces each independent factor to one value
    A = diag(1+rand(nrGens,1));
    b = sign(randn(nrGens,1));
    % instantiate conZonotope with constraints
    conZono = conZonotope(c,G,A,b);
    
    % assert correctness
    if isempty(conZono)
        res = false; break;
    end

    % choose constraints such that conZonotope has to be empty
    % because constraints ||beta|| <= 1 and A*beta = b cannot be fulfilled
    A = diag(0.5*ones(nrGens,1));
    b = sign(randn(nrGens,1));
    % instantiate empty conZonotope
    conZono = conZonotope(c,G,A,b);
    
    % assert correctness
    if ~isempty(conZono)
        res = false; break;
    end
    
end


if res
    disp('testLongDuration_conZonotope_isempty successful');
else
    disp('testLongDuration_conZonotope_isempty failed');
end

%------------- END OF CODE --------------