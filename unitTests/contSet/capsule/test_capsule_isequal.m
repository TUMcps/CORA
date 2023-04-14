function res = test_capsule_isequal
% test_capsule_isequal - unit test function of isequal
%
% Syntax:  
%    res = test_capsule_isequal
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
% Written:      27-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

% empty capsule
C_empty = capsule();

% tolerance
tol = 1e-9;

% define properties
c1 = [2; 0; -1];
c2 = [-1; 1; 0];
g1 = [0.2; -0.7; 0.4];
g2 = [-2; -3; -1];
r1 = 0.2;
r2 = 0.6;
C = capsule(c1,g1,r1);

% test combinations of properties
% ... different center
C_ = capsule(c2,g1,r1);
if isequal(C,C_,tol)
    res = false;
end

% ... different generator
C_ = capsule(c1,g2,r1);
if isequal(C,C_,tol)
    res = false;
end

% ... different radius
C_ = capsule(c1,g1,r2);
if isequal(C,C_,tol)
    res = false;
end

% ... empty capsule
if isequal(C,C_empty)
    res = false;
end

% ... capsule of reduced dimension
C_red = capsule(c1(1:end-1),g1(1:end-1),r1);
if isequal(C,C_red)
    res = false;
end

%------------- END OF CODE --------------