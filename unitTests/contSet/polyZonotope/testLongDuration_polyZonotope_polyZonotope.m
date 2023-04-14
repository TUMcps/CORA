function res = testLongDuration_polyZonotope_polyZonotope
% testLongDuration_polyZonotope_polyZonotope - unit test function of
%    polyZonotope (constructor)
%
% Syntax:  
%    res = testLongDuration_polyZonotope_polyZonotope
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

% Author:       Mark Wetzlinger
% Written:      20-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

tol = 1e-12;

% empty zonotope
pZ = polyZonotope();
res_empty = true;
if ~isempty(pZ)
    res_empty = false;
end

res_rand = true;
nrOfTests = 1000;
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
    Grest = randn(n,nrIndepGens);
    % random exponent matrix (no redudancies) and indentifiers
    counter = 0;
    while true
        counter = counter + 1;
        expMat = randi(7,nrDepFactors,nrDepGens) - 2;
        expMat(expMat < 0) = 0;
        if size(unique(expMat','rows')',2) == size(expMat,2)
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
    if ~all(withinTol(center(pZ),c,tol))
        res_rand = false; break;
    end
    
    % center and dependent generators
    pZ = polyZonotope(c,G);
    if ~all(withinTol(center(pZ),c,tol)) || ~compareMatrices(G,pZ.G,tol)
        res_rand = false; break;
    end
    
    % center and independent generators
    pZ = polyZonotope(c,[],Grest);
    if ~all(withinTol(center(pZ),c,tol)) || ~isempty(pZ.G) ...
            || ~compareMatrices(pZ.Grest,Grest,tol)
        res_rand = false; break;
    end
    
    % center, dependent, and independent generators
    pZ = polyZonotope(c,G,Grest);
    if ~all(withinTol(center(pZ),c,tol)) || ~compareMatrices(pZ.Grest,Grest,tol) ...
            || ~compareMatrices(G,pZ.G,tol)
        res_rand = false; break;
    end
    
    % center, dependent, independent generators, and exponent matrix
    pZ = polyZonotope(c,G,Grest,expMat);
    if ~all(withinTol(center(pZ),c,tol)) || ~compareMatrices(pZ.Grest,Grest,tol) ...
            || ~compareMatrices([G;expMat],[pZ.G;pZ.expMat],tol)
        res_rand = false; break;
    end
    
    % center, dependent generators and exponent matrix
    pZ = polyZonotope(c,G,[],expMat);
    if ~all(withinTol(center(pZ),c,tol)) || ~compareMatrices([G;expMat],[pZ.G;pZ.expMat],tol)
        res_rand = false; break;
    end
    
    % center, dependent generators, exponent matrix, and identifiers
    pZ = polyZonotope(c,G,[],expMat,id);
    if ~all(withinTol(center(pZ),c,tol)) || ~compareMatrices([G;expMat],[pZ.G;pZ.expMat],tol) ...
            || any(abs(pZ.id - id) > tol)
        res_rand = false; break;
    end
    
    % center, independent generators, dependent generators, exponent matrix, and identifiers
    pZ = polyZonotope(c,G,Grest,expMat,id);
    if ~all(withinTol(center(pZ),c,tol)) || ~compareMatrices(pZ.Grest,Grest,tol) ...
            || ~compareMatrices([G;expMat],[pZ.G;pZ.expMat],tol) ...
            || any(abs(pZ.id - id) > tol)
        res_rand = false; break;
    end
    
    % wrong initializations
    % center
    c_mat = randn(n+1);
    c_plus1 = randn(n+1,1);
    randIdx = ceil(n/3);
    c_Inf = c; c_Inf(randIdx) = Inf;
    c_NaN = c; c_NaN(randIdx) = Inf;
    % dependent generators
    G_plus1 = randn(n+1,nrDepGens);
    randLogicals = randn(size(G_plus1)) > 0;
    G_Inf = G_plus1; G_Inf(randLogicals) = Inf;
    G_NaN = G_plus1; G_NaN(randLogicals) = NaN;
    % independent generators
    Grest_plus1 = randn(n+1,nrIndepGens);
    randLogicals = randn(size(Grest_plus1)) > 0;
    Grest_Inf = Grest_plus1; Grest_Inf(randLogicals) = Inf;
    Grest_NaN = Grest_plus1; Grest_NaN(randLogicals) = NaN;
    % exponent matrix
    counter = 0;
    while true
        counter = counter + 1;
        expMat_plus1 = randi(7,nrDepFactors,nrDepGens+1) - 2;
        expMat_plus1(expMat_plus1 < 0) = 0;
        if size(unique(expMat_plus1','rows')',2) == size(expMat_plus1,2)
            break;
        end
        if counter > 10
            continue; % skip current n/nrDepGens/nrDepFactors combination
        end
    end
    % identifiers
    id_plus1 = randperm(nrDepFactors+1)';

    
    % center has Inf/NaN entry
    try
        pZ = polyZonotope(c_Inf); % <- should throw error here
        res_rand = false; break;
    end
    try
        pZ = polyZonotope(c_NaN); % <- should throw error here
        res_rand = false; break;
    end
    
    % center is a matrix
    try
        pZ = polyZonotope(c_mat); % <- should throw error here
        res_rand = false; break;
    end
    
    % empty center
    try
        pZ = polyZonotope([],G); % <- should throw error here
        res_rand = false; break;
    end
    try
        pZ = polyZonotope([],G,Grest); % <- should throw error here
        res_rand = false; break;
    end
    try
        pZ = polyZonotope([],G,Grest,expMat); % <- should throw error here
        res_rand = false; break;
    end
    try
        pZ = polyZonotope([],G,Grest,expMat,id); % <- should throw error here
        res_rand = false; break;
    end
    
    % center and generator matrices of different dimensions
    try
        pZ = polyZonotope(c_plus1,G); % <- should throw error here
        res_rand = false; break;
    end
    try
        pZ = polyZonotope(c_plus1,[],Grest); % <- should throw error here
        res_rand = false; break;
    end
    try
        pZ = polyZonotope(c_plus1,G,Grest); % <- should throw error here
        res_rand = false; break;
    end
    try
        pZ = polyZonotope(c,G_plus1); % <- should throw error here
        res_rand = false; break;
    end
    try
        pZ = polyZonotope(c,[],Grest_plus1); % <- should throw error here
        res_rand = false; break;
    end
    try
        pZ = polyZonotope(c,G,Grest_plus1); % <- should throw error here
        res_rand = false; break;
    end
    
    % dependent generator matrix has Inf/NaN entries
    try
        pZ = polyZonotope(c,G_Inf); % <- should throw error here
        res_rand = false; break;
    end
    try
        pZ = polyZonotope(c,G_NaN); % <- should throw error here
        res_rand = false; break;
    end
    
    % independent generator matrix has Inf/NaN entries
    try
        pZ = polyZonotope(c,[],Grest_Inf); % <- should throw error here
        res_rand = false; break;
    end
    try
        pZ = polyZonotope(c,[],Grest_NaN); % <- should throw error here
        res_rand = false; break;
    end
    
    % exponent matrix and dependent generator matrix do not match
    try
        pZ = polyZonotope(c,G,[],expMat_plus1); % <- should throw error here
        res_rand = false; break;
    end
    try
        pZ = polyZonotope(c,G,[],expMat_plus1,id); % <- should throw error here
        res_rand = false; break;
    end
    try
        pZ = polyZonotope(c,Gplus1,[],expMat); % <- should throw error here
        res_rand = false; break;
    end
    try
        pZ = polyZonotope(c,Gplus1,[],expMat,id); % <- should throw error here
        res_rand = false; break;
    end
    
    % exponent matrix and identifier vector do not match
    try
        pZ = polyZonotope(c,G,[],expMat,id_plus1); % <- should throw error here
        res_rand = false; break;
    end
    
    % too many input arguments
    try
        pZ = polyZonotope(c,G,Grest,expMat,id,id); % <- should throw error here
        res_rand = false; break;
    end 
end


% combine results
res = res_empty && res_rand;

if ~res
    path = pathFailedTests(mfilename());
    save(path,'pZ','c','G','expMat','Grest');
end

%------------- END OF CODE --------------