function res = testLong_zonotope_deleteZeros
% testLong_zonotope_deleteZeros - unit test function of deleteZeros
%    this encompasses checking the function nonzeroFilter
%
% Syntax:  
%    res = testLong_zonotope_deleteZeros
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

% Author:       Matthias Althoff, Mark Wetzlinger
% Written:      26-July-2016
% Last update:  09-August-2020 (MW, enhance randomness)
% Last revision:---

%------------- BEGIN CODE --------------

% assume true
res = true;

% number of tests
nrTests = 1000;

% compare randomly generated zonotopes
for test=1:nrTests

    % random dimension
    n = randi(10);

    % random number of generators
    nrOfGens = randi([2*n,5*n]);
    
    % center
    c = zeros(n,1);
    % generator matrix without any all-zero generators
    Gnozeros = -1+2*rand(n,nrOfGens);
    Znozeros = zonotope(c,Gnozeros);
    
    Zdel1 = deleteZeros(Znozeros);
    
    % since no zero generators, results has to be the same as before
    res1 = compareMatrices(c,center(Zdel1)) ...
        && compareMatrices(generators(Zdel1),Gnozeros);
    
    % insert zero generators in matrix at random position
    Gzeros = [Gnozeros,zeros(n,randi([5,25],1,1))];
    Gzeros = Gzeros(:,randperm(size(Gzeros,2)));
    Zzeros = zonotope(c,Gzeros);
    
    % delete zero generators
    Zdel2 = deleteZeros(Zzeros);
    
    % result has to be the same as original zonotope
    res2 = isequal(Znozeros,Zdel2);

    % add results
    if ~(res1 && res2)
        res = false; return
    end
    
end

%------------- END OF CODE --------------
