function res = testLongDuration_conZonotope_randPoint
% testLongDuration_conZonotope_randPoint - unit test function of randPoint
%
% Syntax:  
%    res = testLongDuration_conZonotope_randPoint
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
tol = 1e-9;

% check empty conZonotope object -> error
conZono = conZonotope();
try 
    p = randPoint(conZono);
catch ME
    if ~strcmp(ME.identifier,'CORA:emptySet')
        res = false;
    end
end


% number of tests
nrOfTests = 100;

for i=1:nrOfTests
    % random dimension
    n = randi([2,8]); % small because of containment checks
    
    % random center
    c = randn(n,1);
    % random generator matrix
    G = randn(n,n+randi(5));
    nrGens = size(G,2);
    
    % instantiate conZonotope without constraints
    conZono = conZonotope(c,G);
    
    % compute random points
    nrPts = 10;
    pNormal = randPoint(conZono,nrPts,'standard');
    pExtreme = randPoint(conZono,nrPts,'extreme');
    
    % check for containment in zonotope
    Z = zonotope(c,G);
    res = in(Z,pNormal,tol);
    
    % random constraints so that conZonotope represents just a point
    % as A being diagional forces each independent factor to one value
    A = diag(1+rand(nrGens,1));
    b = sign(randn(nrGens,1));
    % instantiate conZonotope with constraints
    conZono = conZonotope(c,G,A,b);
    
    % compute single possible point -> bug in randPoint
%     p = randPoint(conZono,1,'standard');
    
    % compute point analytically
    p_analytical = zeros(n,1);
    conZonoG = conZono.Z(:,2:end);
    beta = zeros(nrGens,1);
    for k=1:nrGens
        % constraints yield only one possible value for beta_k
        beta(k) = conZono.b(k) / conZono.A(k,k);
        
        % check if beta_k in [-1,1]
        if abs(beta(k)) > 1
            error("Check test function!");
        end
        
        % contribution of k-th generator
        kthGen = beta(k) * conZonoG(:,k);
        % add to point
        p_analytical = p_analytical + kthGen;
    end
    p = p_analytical;
    
    % compare results
    if any(abs(p - p_analytical) > tol)
        res = false; break;
    end
    
end


if res
    disp('testLongDuration_conZonotope_randPoint successful');
else
    disp('testLongDuration_conZonotope_randPoint failed');
end

%------------- END OF CODE --------------