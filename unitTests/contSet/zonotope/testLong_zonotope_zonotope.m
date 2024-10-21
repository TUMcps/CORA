function res = testLong_zonotope_zonotope
% testLong_zonotope_zonotope - unit test function of zonotope (constructor)
%
% Syntax:
%    res = testLong_zonotope_zonotope
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
% Written:       20-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% number of tests
nrOfTests = 1000;

for i=1:nrOfTests

    % random dimension
    n = randi(25);
    % random number of generators
    nrGens = randi(25);
    
    % random center, random generator matrix
    c = randn(n,1);
    G = randn(n,nrGens);
    Gempty = [];
    Zmat = [c,G];
    
    % admissible initializations
    Z = zonotope(c,G);
    assertLoop(compareMatrices(center(Z),c),i)
    assertLoop(compareMatrices(generators(Z),G),i)
    
    Z = zonotope(c,Gempty);
    assertLoop(compareMatrices(center(Z),c),i)
    assertLoop(isempty(generators(Z)),i)
    
    Z = zonotope(Zmat);
    assertLoop(compareMatrices(center(Z),Zmat(:,1)),i)
    assertLoop(compareMatrices(generators(Z),Zmat(:,2:end)),i)
    
    % wrong initializations
    c_plus1 = randn(n+1,1);
    G_plus1 = randn(n+1,nrGens);
    randIdx = ceil(n/3);
    c_NaN = c; c_NaN(randIdx) = NaN;
    G_NaN = G; 
    while ~any(isnan(G_NaN))
        randLogicals = randn(size(G)) > 0;
        G_NaN(randLogicals) = NaN;
    end
    
    % center and generator matrix do not match
    assertThrowsAs(@zonotope,'CORA:wrongInputInConstructor',c_plus1,G);
    assertThrowsAs(@zonotope,'CORA:wrongInputInConstructor',c,G_plus1);
    
    % center is empty
    assertThrowsAs(@zonotope,'CORA:wrongInputInConstructor',[],G);
    
    % center has NaN entry
    assertThrowsAs(@zonotope,'CORA:wrongValue',c_NaN,G);
    
    % generator matrix has NaN entries
    assertThrowsAs(@zonotope,'CORA:wrongValue',c,G_NaN);
    
    % too many input arguments
    assertThrowsAs(@zonotope,'CORA:numInputArgsConstructor',c,G,G);
    
end

% ------------------------------ END OF CODE ------------------------------
