function res = testLongDuration_zonotope_zonotope
% testLongDuration_zonotope_zonotope - unit test function of zonotope (constructor)
%
% Syntax:  
%    res = testLongDuration_zonotope_zonotope
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
res = true;

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
    if any(abs(center(Z) - c) > tol) || any(any(abs(generators(Z) - G) > tol))
        res = false; break;
    end
    
    Z = zonotope(c,Gempty);
    if any(abs(center(Z) - c) > tol) || ~isempty(generators(Z))
        res = false; break;
    end
    
    Z = zonotope(Zmat);
    if any(abs(center(Z) - Zmat(:,1)) > tol) ...
            || any(any(abs(generators(Z) - Zmat(:,2:end)) > tol))
        res = false; break;
    end
    
    % wrong initializations
    c_plus1 = randn(n+1,1);
    G_plus1 = randn(n+1,nrGens);
%     randIdx = ceil(n/3);
%     randLogicals = randn(size(G)) > 0;
%     c_Inf = c; c_Inf(randIdx) = Inf;
%     c_NaN = c; c_NaN(randIdx) = Inf;
%     G_Inf = G; G_Inf(randLogicals) = Inf;
%     G_NaN = G; G_NaN(randLogicals) = NaN;
    
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
    
    
    % center has Inf/NaN entry
%     try
%         Z = zonotope(c_Inf,G); % <- should throw error here
%         res = false; break;
%     end
%     try
%         Z = zonotope(c_NaN,G); % <- should throw error here
%         res = false; break;
%     end
    
    % generator matrix has Inf/NaN entries
%     try
%         Z = zonotope(c,G_Inf); % <- should throw error here
%         res = false; break;
%     end
%     try
%         Z = zonotope(c,G_NaN); % <- should throw error here
%         res = false; break;
%     end
    
    % too many input arguments
    try
        Z = zonotope(c,G,G); % <- should throw error here
        res = false; break;
    end 
end


if res
    disp('testLongDuration_zonotope successful');
else
    disp('testLongDuration_zonotope failed');
end

%------------- END OF CODE --------------