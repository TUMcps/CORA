function res = test_capsule_capsule
% test_capsule_capsule - unit test function of capsule (constructor)
%
% Syntax:  
%    res = test_capsule_capsule
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

% assume true
res = true;

% empty capsule
C_empty = capsule();
if ~isempty(center(C_empty))
    res = false;
end


% random center, generator, and radius
c = [2; 0];
g = [1; -1];
r = 0.5;

% admissible initializations
% only center
C = capsule(c);
if ~compareMatrices(C.c,c)
    res = false;
end

% center and generator
C = capsule(c,g);
if ~compareMatrices(C.c,c) || ~compareMatrices(C.g,g)
    res = false;
end

% center and radius:
C = capsule(c,r);
if ~compareMatrices(C.c,c) || ~withinTol(C.r,r)
	res = false;
end

% center, generator, and radius
C = capsule(c,g,r);
if ~compareMatrices(C.c,c) || ~compareMatrices(C.g,g) || ~withinTol(C.r,r)
    res = false;
end

% wrong initializations
cn_1 = [3; 2; -1];
gn_1 = [-2; 1; 1];
rneg = -0.2;
rvec = [1; 4];

% mismatch between center and generator
try
    C = capsule(c,gn_1); % <- should throw error here
    res = false;
end
try
    C = capsule(cn_1,g); % <- should throw error here
    res = false;
end

% negative radius
try
    C = capsule(c,g,rneg); % <- should throw error here
    res = false;
end

% radius as a vector
try
    C = capsule(c,g,rvec); % <- should throw error here
    res = false;
end

% too many input arguments
try
    C = capsule(c,g,r,r); % <- should throw error here
    res = false;
end

%------------- END OF CODE --------------