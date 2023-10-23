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
% See also: -

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
    if ~compareMatrices(center(Z),c) || ~compareMatrices(generators(Z),G)
        res = false; break;
    end
    
    Z = zonotope(c,Gempty);
    if ~compareMatrices(center(Z),c) || ~isempty(generators(Z))
        res = false; break;
    end
    
    Z = zonotope(Zmat);
    if ~compareMatrices(center(Z),Zmat(:,1)) ...
            || ~compareMatrices(generators(Z),Zmat(:,2:end))
        res = false; break;
    end
    
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
    try
        Z = zonotope(c_plus1,G); % <- should throw error here
        res = false; break;
    end
    try
        Z = zonotope(c,G_plus1); % <- should throw error here
        res = false; break;
    end
    
    % center is empty
    try
        Z = zonotope([],G); % <- should throw error here
        res = false; break;
    end
    
    
    % center has NaN entry
    try
        Z = zonotope(c_NaN,G); % <- should throw error here
        res = false; break;
    end
    
    % generator matrix has NaN entries
    try
        Z = zonotope(c,G_NaN); % <- should throw error here
        res = false; break;
    end
    
    % too many input arguments
    try
        Z = zonotope(c,G,G); % <- should throw error here
        res = false; break;
    end 
end

% ------------------------------ END OF CODE ------------------------------
