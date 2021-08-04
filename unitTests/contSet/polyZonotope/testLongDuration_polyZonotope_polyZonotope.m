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
    if any(abs(center(pZ) - c) > tol)
        res_rand = false; break;
    end
    
    % center and dependent generators
    pZ = polyZonotope(c,G);
    if any(abs(center(pZ) - c) > tol) || any(any(abs(pZ.G - G) > tol))
        res_rand = false; break;
    end
    
    % center and independent generators
    pZ = polyZonotope(c,[],Grest);
    if any(abs(center(pZ) - c) > tol) || ~isempty(pZ.G) ...
            || any(any(abs(pZ.Grest - Grest) > tol))
        res_rand = false; break;
    end
    
    % center, dependent, and independent generators
    pZ = polyZonotope(c,G,Grest);
    if any(abs(center(pZ) - c) > tol) || any(any(abs(pZ.G - G) > tol)) ...
            || any(any(abs(pZ.Grest - Grest) > tol))
        res_rand = false; break;
    end
    
    % center, dependent, independent generators, and exponent matrix
    pZ = polyZonotope(c,G,Grest,expMat);
    if any(abs(center(pZ) - c) > tol) || any(any(abs(pZ.Grest - Grest) > tol)) ...
            || ~sameGexpMat(G,expMat,pZ.G,pZ.expMat,tol)
        res_rand = false; break;
    end
    
    % center, dependent generators and exponent matrix
    pZ = polyZonotope(c,G,[],expMat);
    if any(abs(center(pZ) - c) > tol) || ~sameGexpMat(G,expMat,pZ.G,pZ.expMat,tol)
        res_rand = false; break;
    end
    
    % center, dependent generators, exponent matrix, and identifiers
    pZ = polyZonotope(c,G,[],expMat,id);
    if any(abs(center(pZ) - c) > tol) || ~sameGexpMat(G,expMat,pZ.G,pZ.expMat,tol) ...
            || any(abs(pZ.id - id) > tol)
        res_rand = false; break;
    end
    
    % center, independent generators, dependent generators, exponent matrix, and identifiers
    pZ = polyZonotope(c,G,Grest,expMat,id);
    if any(abs(center(pZ) - c) > tol) || any(any(abs(pZ.Grest - Grest) > tol)) ...
            || ~sameGexpMat(G,expMat,pZ.G,pZ.expMat,tol) ...
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

if res
    disp('testLongDuration_polyZonotope successful');
else
    disp('testLongDuration_polyZonotope failed');
end

end


% Auxiliary Function ------------------------------------------------------
function same = sameGexpMat(G,expMat,G_new,expMat_new,tol)
% assumption: no redundancies removal in expMat -> expMat_new

if any(size(G) ~= size(G_new)) || any(size(expMat) ~= size(expMat_new))
    same = false; return;
end

% quick check (same order as before)
if ~(any(any(abs(G_new - G) > tol)) || any(any(abs(expMat_new - expMat) > tol)))
    same = true; return;
end

% init output argument
same = true;
% go through generators and exponent columns one by one
nrGens = size(G,2);
for k=1:nrGens
    % select first generator
    gk = G(:,1);
    
    % search for gk in G_new
    idx = find(all(abs(G_new - gk) < tol,1),1,'first');
    % generator not found
    if isempty(idx)
        same = false; break;
    end
    
    % check equality of expMat columns
    if any(abs(expMat(:,1) - expMat_new(:,idx)) > tol)
        same = false; break;
    end
    
    % remove columns in Gs and expMats
    G(:,1) = []; expMat(:,1) = [];
    G_new(:,idx) = []; expMat_new(:,idx) = [];
end

end

%------------- END OF CODE --------------