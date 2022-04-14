function res = test_capsule_randPoint
% test_capsule_randPoint - unit test function of randPoint
%
% Syntax:  
%    res = test_capsule_randPoint
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required:
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Mark Wetzlinger
% Written:      27-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% empty case
C = capsule();
res(1) = true;
try
    p = randPoint(C); % <- should throw error here
    res(1) = false;
catch ME
    if ~strcmp(ME.identifier,'CORA:emptySet')
        res(1) = false;
    end
end

% degenerate capsule
c = [1; -1; 2];
g = [0; 0; 0];
r = 0;
C = capsule(c,g,r);

% rand point only center
p = randPoint(C);

% check
res(2) = all(abs(p - c) < 5*eps);

% instantiate full-dimensional capsule
c = [3; 0; 0];
g = [-1; 3; 2];
r = 1;
C = capsule(c,g,r);

% generate random points
numPoints = 10;
p = zeros(dim(C),numPoints);
for i=1:numPoints
    p(:,i) = randPoint(C);
end

% check if all random points inside capsule
res(3) = in(C,p);

% combine results
res = all(res);

if res
    disp('test_randPoint successful');
else
    disp('test_randPoint failed');
end

%------------- END OF CODE --------------