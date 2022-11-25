function res = testLongDuration_conZonotope_project
% testLongDuration_conZonotope_project - unit test function of project
%
% Syntax:  
%    res = testLongDuration_conZonotope_project
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

% check empty conZonotope object: projection should return error
conZono = conZonotope();
try 
    project(conZono,[1,2])
    res = false;
catch ME
    if ~strcmp(ME.identifier,'CORA:emptySet')
        res = false;
    end
end


% number of tests
nrOfTests = 1000;

for i=1:nrOfTests
    % random dimension
    n = randi([2,50]);
    % projection to random subspace
    selectedDims = randperm(n);
    selectedDims = selectedDims(1:floor(n/2));
    selectedDims = sort(selectedDims);
    projDim = 1:n;
    projDim = projDim(selectedDims);
    
    % random center
    c = randn(n,1);
    % random generator matrix
    G = randn(n,randi(10));
    nrGens = size(G,2);
    
    % instantiate conZonotope without constraints
    conZono = conZonotope(c,G);
    
    % projected conZonotope
    conZonoProj = project(conZono,projDim);
    
    % instantiate projection
    conZonoProj_true = conZonotope(c(selectedDims),G(selectedDims,:));
       
    % compare results
    if any(any(abs(conZonoProj.Z - conZonoProj_true.Z))) > tol || ...
            any(any(abs(conZonoProj.A - conZonoProj_true.A))) > tol || ...
            any(any(abs(conZonoProj.b - conZonoProj_true.b))) > tol
        res = false; break;
    end
    
    % random constraints so that conZonotope represents just a point
    % as A being diagional forces each independent factor to one value
    A = diag(1+rand(nrGens,1));
    b = sign(randn(nrGens,1));
    % instantiate conZonotope with constraints
    conZono = conZonotope(c,G,A,b);
    
    % project to subspace
    conZonoProj = project(conZono,projDim);
    
    % instantiate projection
    conZonoProj_true = conZonotope(c(selectedDims),G(selectedDims,:),A,b);
    
    % compare results
    if any(any(abs(conZonoProj.Z - conZonoProj_true.Z))) > tol || ...
            any(any(abs(conZonoProj.A - conZonoProj_true.A))) > tol || ...
            any(any(abs(conZonoProj.b - conZonoProj_true.b))) > tol
        res = false; break;
    end
    
end


if res
    disp('testLongDuration_conZonotope_project successful');
else
    disp('testLongDuration_conZonotope_project failed');
end

%------------- END OF CODE --------------