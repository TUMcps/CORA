function res = test_zonotope_containsPoint
% test_zonotope_containsPoint - unit test function of containsPoint
%
% Syntax:  
%    res = test_zonotope_containsPoint
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
% Written:      26-August-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% create a zonotope
dim = 2;
gens = 5;
% rand gives value of [0,1]
Z = zonotope([zeros(dim,1),rand(dim,gens)]);

% Single points -----------------------------------------------------------
% create a point inside: center of zonotope
p_inside = center(Z);
% ...and outside: add all generators (with max of rand -> 1)
p_outside = gens*ones(dim,1);

% check if correct results for containment
res_inside  = containsPoint(Z, p_inside);
res_outside = containsPoint(Z, p_outside);

% point on the border should also count as outside (?)
% p_border = sum(generators(Z),2);
% res_border = containsPoint(Z, p_border);

res_single = res_inside && ~res_outside;% && ~res_border;
% -------------------------------------------------------------------------

% Array of points ---------------------------------------------------------
% all outside...
num = 10;
p_array = gens*(ones(dim,num)+rand(dim,num));
res_array = any(containsPoint(Z,p_array));
% -------------------------------------------------------------------------

% add results
res = res_single && ~res_array;

if res
    disp('test_containsPoint successful');
else
    disp('test_containsPoint failed');
end

%------------- END OF CODE --------------