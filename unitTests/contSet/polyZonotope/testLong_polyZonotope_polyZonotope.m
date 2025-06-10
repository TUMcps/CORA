function res = testLong_polyZonotope_polyZonotope
% testLong_polyZonotope_polyZonotope - unit test function of
%    polyZonotope (constructor)
%
% Syntax:
%    res = testLong_polyZonotope_polyZonotope
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
 
% tolerance
tol = 1e-12;

nrOfTests = 50;
for i=1:nrOfTests

    % random dimension
    n = randi(10);
    % random number of generators
    nrIndepGens = randi(25);
    nrDepGens = randi(25);
    % random number of dependent factors
    nrDepFactors = randi([nrDepGens,nrDepGens+10]);
    
    % random center, random generator matrices
    c = randn(n,1);
    G = randn(n,nrDepGens);
    GI = randn(n,nrIndepGens);
    % random exponent matrix (no redudancies) and indentifiers
    counter = 0;
    while true
        counter = counter + 1;
        E = randi(7,nrDepFactors,nrDepGens) - 2;
        E(E < 0) = 0;
        if size(unique(E','rows')',2) == size(E,2)
            break;
        end
        if counter > 10
            continue; % skip current n/nrDepGens/nrDepFactors combination
        end
    end
    id = randperm(nrDepFactors)';
    
    % admissible initializations
    % only center
    pZ = polyZonotope(c);
    assert(all(withinTol(center(pZ),c,tol)));
    
    % center and dependent generators
    pZ = polyZonotope(c,G);
    assert(all(withinTol(center(pZ),c,tol)));
    assert(compareMatrices(G,pZ.G,tol));
    
    % center and independent generators
    pZ = polyZonotope(c,[],GI);
    assert(all(withinTol(center(pZ),c,tol)));
    assert(isempty(pZ.G));
    assert(compareMatrices(pZ.GI,GI,tol));
    
    % center, dependent, and independent generators
    pZ = polyZonotope(c,G,GI);
    assert(all(withinTol(center(pZ),c,tol)));
    assert(compareMatrices(pZ.GI,GI,tol));
    assert(compareMatrices(G,pZ.G,tol));
    
    % center, dependent, independent generators, and exponent matrix
    pZ = polyZonotope(c,G,GI,E);
    assert(all(withinTol(center(pZ),c,tol)));
    assert(compareMatrices(pZ.GI,GI,tol));
    assert(compareMatrices([G;E],[pZ.G;pZ.E],tol));
    
    % center, dependent generators and exponent matrix
    pZ = polyZonotope(c,G,[],E);
    assert(all(withinTol(center(pZ),c,tol)))
    assert(compareMatrices([G;E],[pZ.G;pZ.E],tol))
    
    % center, dependent generators, exponent matrix, and identifiers
    pZ = polyZonotope(c,G,[],E,id);
    assert(all(withinTol(center(pZ),c,tol)));
    assert(compareMatrices([G;E],[pZ.G;pZ.E],tol));
    assert(all(abs(pZ.id - id) <= tol));
    
    % center, independent generators, dependent generators, exponent matrix, and identifiers
    pZ = polyZonotope(c,G,GI,E,id);
    assert(all(withinTol(center(pZ),c,tol)));
    assert(compareMatrices(pZ.GI,GI,tol));
    assert(compareMatrices([G;E],[pZ.G;pZ.E],tol));
    assert(all(abs(pZ.id - id) <= tol));
    
    % wrong initializations
    % center
    c_mat = randn(n+1);
    c_plus1 = randn(n+1,1);
    randIdx = ceil(n/3);
    c_Inf = c; c_Inf(randIdx) = Inf;
    c_NaN = c; c_NaN(randIdx) = Inf;
    % dependent generators
    G_plus1 = randn(n+1,nrDepGens);
    G_Inf = G; G_Inf(randi(numel(G_Inf))) = Inf;
    G_NaN = G; G_NaN(randi(numel(G_NaN))) = NaN;
    % independent generators
    GI_plus1 = randn(n+1,nrIndepGens);
    GI_Inf = G; GI_Inf(randi(numel(GI_Inf))) = Inf;
    GI_NaN = G; GI_NaN(randi(numel(GI_NaN))) = NaN;
    % exponent matrix
    counter = 0;
    while true
        counter = counter + 1;
        E_plus1 = randi(7,nrDepFactors,nrDepGens+1) - 2;
        E_plus1(E_plus1 < 0) = 0;
        if size(unique(E_plus1','rows')',2) == size(E_plus1,2)
            break;
        end
        if counter > 10
            continue; % skip current n/nrDepGens/nrDepFactors combination
        end
    end
    % identifiers
    id_plus1 = randperm(nrDepFactors+1)';

    
    % center has Inf/NaN entry
    assertThrowsAs(@polyZonotope,'CORA:wrongValue',c_Inf);
    assertThrowsAs(@polyZonotope,'CORA:wrongValue',c_NaN);

    % center is a matrix
    assertThrowsAs(@polyZonotope,'CORA:wrongInputInConstructor',c_mat);

    % empty center
    assertThrowsAs(@polyZonotope,'CORA:wrongInputInConstructor',[],G);
    assertThrowsAs(@polyZonotope,'CORA:wrongInputInConstructor',[],G,GI);
    assertThrowsAs(@polyZonotope,'CORA:wrongInputInConstructor',[],G,GI,E);
    assertThrowsAs(@polyZonotope,'CORA:wrongInputInConstructor',[],G,GI,E,id);

    % center and generator matrices of different dimensions
    assertThrowsAs(@polyZonotope,'CORA:wrongInputInConstructor',c_plus1,G);
    assertThrowsAs(@polyZonotope,'CORA:wrongInputInConstructor',c_plus1,[],GI);
    assertThrowsAs(@polyZonotope,'CORA:wrongInputInConstructor',c_plus1,G,GI);
    assertThrowsAs(@polyZonotope,'CORA:wrongInputInConstructor',c,G_plus1);
    assertThrowsAs(@polyZonotope,'CORA:wrongInputInConstructor',c,[],GI_plus1);
    assertThrowsAs(@polyZonotope,'CORA:wrongInputInConstructor',c,G,GI_plus1);
    
    % dependent generator matrix has Inf/NaN entries
    assertThrowsAs(@polyZonotope,'CORA:wrongValue',c,G_Inf);
    assertThrowsAs(@polyZonotope,'CORA:wrongValue',c,G_NaN);

    % independent generator matrix has Inf/NaN entries
    assertThrowsAs(@polyZonotope,'CORA:wrongValue',c,[],GI_Inf);
    assertThrowsAs(@polyZonotope,'CORA:wrongValue',c,[],GI_NaN);
    
    % exponent matrix and dependent generator matrix do not match
    assertThrowsAs(@polyZonotope,'CORA:wrongInputInConstructor',c,G,[],E_plus1);
    assertThrowsAs(@polyZonotope,'CORA:wrongInputInConstructor',c,G,[],E_plus1,id);
    assertThrowsAs(@polyZonotope,'CORA:wrongInputInConstructor',c,G_plus1,[],E);
    assertThrowsAs(@polyZonotope,'CORA:wrongInputInConstructor',c,G_plus1,[],E,id);
    
    % exponent matrix and identifier vector do not match
    assertThrowsAs(@polyZonotope,'CORA:wrongInputInConstructor',c,G,[],E,id_plus1);
    
    % too many input arguments
    assertThrowsAs(@polyZonotope,'CORA:numInputArgsConstructor',c,G,GI,E,id,id);
     
end

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
