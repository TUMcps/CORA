function res = testLong_conZonotope_deleteZeros
% testLong_conZonotope_deleteZeros - unit test function for the deletion of
%    redundant generators or constraints
%
% Syntax:
%    res = testLong_conZonotope_deleteZeros
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
% Written:       28-March-2022
% Last update:   ---
% Last revision: 09-January-2024 (MW, rename testLong_*)

% ------------------------------ BEGIN CODE -------------------------------

res = true;

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
    nrGensRemoved = length(idx);
    Gzeros(:,idx) = 0;

    % random constraints
    nrConstr = randi([1,10]);
    A = randn(nrConstr,nrGens);
    b = randn(nrConstr,1);

    % set random indices to zero
    idx = unique(randi([1,nrConstr],ceil(nrConstr/4),1));
    nrConstrRemoved = length(idx);
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
    cZ_noG_ = deleteZeros(cZ_noG);
    cZ_noAb_ = deleteZeros(cZ_noAb);
    cZ_Gzeros_noAb_ = deleteZeros(cZ_Gzeros_noAb);
    cZ_ = deleteZeros(cZ);
    cZ_Abzeros_ = deleteZeros(cZ_Abzeros);
    cZ_GAbzeros_ = deleteZeros(cZ_GAbzeros);

    % checks
    % no generator matrix (no A,b) -> no change
    if size(cZ_noG_.G,2) ~= 0
        res = false; break;
    end
    % full G, no A,b -> no change
    if size(cZ_noAb_.G,2) ~= nrGens
        res = false; break;
    end
    % full G,A,b -> no change
    if size(cZ_.G,2) ~= nrGens
        res = false; break;
    end
    if size(cZ_.A,1) ~= nrConstr
        res = false; break;
    end
    % G with zeros, no A,b -> less gens
    if size(cZ_Gzeros_noAb_.G,2) ~= nrGens - nrGensRemoved
        res = false; break;
    end
    % full G, A,b with zeros -> less constr
    if size(cZ_Abzeros_.G,2) ~= nrGens
        res = false; break;
    end
    if size(cZ_Abzeros_.A,1) ~= nrConstr - nrConstrRemoved
        res = false; break;
    end
    % G,A,b with zeros -> less gens, less constr
    if size(cZ_GAbzeros_.A,1) ~= nrConstr - nrConstrRemoved
        res = false; break;
    end

end


% ------------------------------ END OF CODE ------------------------------
