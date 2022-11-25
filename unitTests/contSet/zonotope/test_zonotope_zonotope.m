function res = test_zonotope_zonotope
% test_zonotope_zonotope - unit test function of zonotope (constructor)
%
% Syntax:  
%    res = test_zonotope_zonotope
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
% Written:      27-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

tol = 1e-12;

% empty zonotope
Z = zonotope();
res = isempty(Z);


% random center, random generator matrix
c = [3; 3; 2];
G = [2 -4 -6 3 5; 1 -7 3 -5 2; 0 4 -7 3 2];
Gempty = [];
Zmat = [c,G];

% admissible initializations
Z = zonotope(c,G);
if any(abs(center(Z) - c) > tol) || any(any(abs(generators(Z) - G) > tol))
    res = false;
end

Z = zonotope(c,Gempty);
if any(abs(center(Z) - c) > tol) || ~isempty(generators(Z))
    res = false;
end

Z = zonotope(Zmat);
if any(abs(center(Z) - Zmat(:,1)) > tol) ...
        || any(any(abs(generators(Z) - Zmat(:,2:end)) > tol))
    res = false;
end

% wrong initializations
c_plus1 = [4; 6; -2; 3];
G_plus1 = [2 -4 -6 3 5; 1 -7 3 -5 2; 0 4 -7 3 2; 2 0 5 -4 2];
% randIdx = ceil(n/3);
% randLogicals = randn(size(G)) > 0;
% c_Inf = c; c_Inf(randIdx) = Inf;
% c_NaN = c; c_NaN(randIdx) = Inf;
% G_Inf = G; G_Inf(randLogicals) = Inf;
% G_NaN = G; G_NaN(randLogicals) = NaN;

% center and generator matrix do not match
try
    Z = zonotope(c_plus1,G); % <- should throw error here
    res = false;
end
try
    Z = zonotope(c,G_plus1); % <- should throw error here
    res = false;
end

% center is empty
try
    Z = zonotope([],G); % <- should throw error here
    res = false;
end


% center has Inf/NaN entry
% try
%     Z = zonotope(c_Inf,G); % <- should throw error here
%     res = false;
% end
% try
%     Z = zonotope(c_NaN,G); % <- should throw error here
%     res = false;
% end

% generator matrix has Inf/NaN entries
% try
%     Z = zonotope(c,G_Inf); % <- should throw error here
%     res = false;
% end
% try
%     Z = zonotope(c,G_NaN); % <- should throw error here
%     res = false;
% end

% too many input arguments
try
    Z = zonotope(c,G,G); % <- should throw error here
    res = false;
end 



if res
    disp('test_zonotope successful');
else
    disp('test_zonotope failed');
end

%------------- END OF CODE --------------