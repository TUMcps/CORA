function res = test_conZonotope_generators
% test_conZonotope_generators - unit test function for reading out the
%    generator matrix of a conZonotope object
%
% Syntax:  
%    res = test_conZonotope_generators
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
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Author:       Mark Wetzlinger
% Written:      28-March-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

% check empty conZonotope
cZ_empty = conZonotope();
cZ_empty_ = generators(cZ_empty);
if ~isempty(cZ_empty_)
    res = false; return;
end


% number of tests
nrOfTests = 100;

for i=1:nrOfTests
    % random dimension
    n = randi([1,10]);
    % random center
    c = randn(n,1);
    % random generator matrix (at least n generators)
    G = randn(n,n+randi(10));
    nrGens = size(G,2);

    % set random indices to zero
    Gzeros = G;
    idx = unique(randi([1,nrGens],ceil(nrGens/4),1));
    Gzeros(:,idx) = 0;

    % random constraints
    nrConstr = randi([1,10]);
    A = randn(nrConstr,nrGens);
    b = randn(nrConstr,1);

    % set random indices to zero
    idx = unique(randi([1,nrConstr],ceil(nrConstr/4),1));
    Azeros = A;
    Azeros(idx,:) = 0;
    bzeros = b;
    bzeros(idx) = 0;

    % instantiate conZonotope without generator matrix
    cZ_noG = conZonotope(c);

    % instantiate conZonotope without constraints
    cZ_noAb = conZonotope(c,G);
    cZ_Gzeros_noAb = conZonotope(c,Gzeros);

    % instantiate conZonotope with constraints
    cZ = conZonotope(c,G,A,b);
    cZ_Abzeros = conZonotope(c,G,Azeros,bzeros);
    cZ_GAbzeros = conZonotope(c,Gzeros,Azeros,bzeros);

    % remove zero-length generators
    cZ_noG_ = generators(cZ_noG);
    cZ_noAb_ = generators(cZ_noAb);
    cZ_Gzeros_noAb_ = generators(cZ_Gzeros_noAb);
    cZ_ = generators(cZ);
    cZ_Abzeros_ = generators(cZ_Abzeros);
    cZ_GAbzeros_ = generators(cZ_GAbzeros);

    % checks
    % no generator matrix (no A,b)
    if size(cZ_noG_,2) ~= 0
        res = false; break;
    end
    % full G, no A,b
    if size(cZ_noAb_,2) ~= nrGens || ~compareMatrices(G,cZ_noAb_)
        res = false; break;
    end
    % full G,A,b
    if size(cZ_,2) ~= nrGens || ~compareMatrices(G,cZ_)
        res = false; break;
    end
    % G with zeros, no A,b
    if size(cZ_Gzeros_noAb_,2) ~= nrGens || ~compareMatrices(Gzeros,cZ_Gzeros_noAb_)
        res = false; break;
    end
    % full G, A,b with zeros
    if size(cZ_Abzeros_,2) ~= nrGens || ~compareMatrices(G,cZ_Abzeros_)
        res = false; break;
    end
    % G,A,b with zeros
    if size(cZ_GAbzeros_,2) ~= nrGens || ~compareMatrices(Gzeros,cZ_GAbzeros_)
        res = false; break;
    end

end


%------------- END OF CODE --------------