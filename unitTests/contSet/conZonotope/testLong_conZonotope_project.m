function res = testLong_conZonotope_project
% testLong_conZonotope_project - unit test function of project
%
% Syntax:
%    res = testLong_conZonotope_project
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
% Written:       14-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% check empty conZonotope object: projection should return error
cZ = conZonotope.empty(1);
try 
    project(cZ,[1,2])
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
    if rand > 0.5
        % use logical indexing
        projDim = false(n,1);
        projDim(selectedDims) = true;
    else
        % use standard indexing
        projDim = 1:n;
        projDim = projDim(selectedDims);
    end
    
    % random center
    c = randn(n,1);
    % random generator matrix
    G = randn(n,randi(10));
    nrGens = size(G,2);
    
    % instantiate conZonotope without constraints
    cZ = conZonotope(c,G);
    
    % projected conZonotope
    cZ_proj = project(cZ,projDim);
    
    % instantiate projection
    cZ_proj_true = conZonotope(c(selectedDims),G(selectedDims,:));
       
    % compare results
    if ~isequal(cZ_proj,cZ_proj_true)
        res = false; break;
    end
    
    % random constraints so that conZonotope represents just a point
    % as A being diagional forces each independent factor to one value
    A = diag(1+rand(nrGens,1));
    b = sign(randn(nrGens,1));
    % instantiate conZonotope with constraints
    cZ = conZonotope(c,G,A,b);
    
    % project to subspace
    cZ_proj = project(cZ,projDim);
    
    % instantiate projection
    cZ_proj_true = conZonotope(c(selectedDims),G(selectedDims,:),A,b);
    
    % compare results
    if ~isequal(cZ_proj,cZ_proj_true)
        res = false; break;
    end
    
end

% ------------------------------ END OF CODE ------------------------------
