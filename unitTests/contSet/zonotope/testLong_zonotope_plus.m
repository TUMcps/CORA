function res = testLong_zonotope_plus
% testLong_zonotope_plus - unit test function of plus
%
% Syntax:
%    res = testLong_zonotope_plus
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

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       26-July-2016
% Last update:   09-August-2020
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% number of tests
nrTests = 1000;

for i=1:nrTests

    % random dimension
    n = randi(10);

    % create a random zonotope
    nrOfGens = randi([10,25],1,1);
    c1 = -1+2*rand(n,1);
    G1 = -1+2*rand(n,nrOfGens);
    Z1 = zonotope(c1,G1);
    c2 = -1+2*rand(n,1);
    G2 = -1+2*rand(n,nrOfGens);
    Z2 = zonotope(c2,G2);
    
    % add zonotopes
    Zres = Z1 + Z2;

    % check center
    if ~compareMatrices(c1+c2,center(Zres))
        res = false; return
    end
    
    % check generators
    Gres = generators(Zres);
    G1inGres = nnz(all(ismember(Gres,G1),1));
    G2inGres = nnz(all(ismember(Gres,G2),1));
    
    if size(Gres,2) ~= G1inGres+G2inGres
        res = false; return
    end
end

% ------------------------------ END OF CODE ------------------------------
