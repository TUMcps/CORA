function res = test_capsule_isZero
% test_capsule_isZero - unit test function of isZero
%
% Syntax:  
%    res = test_capsule_isZero
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

% Author:       Mark Wetzlinger
% Written:      16-March-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% empty case
res(1) = ~isZero(capsule());

% check different cases...

% true isZero cases
C = capsule(zeros(3,1));
res(end+1) = isZero(C);
C = capsule(zeros(3,1),zeros(3,1),0);
res(end+1) = isZero(C);

% shifted center
C = capsule(ones(3,1),zeros(3,1),0);
res(end+1) = ~isZero(C);

% including generator, no radius
C = capsule(zeros(2,1),ones(2,1),0);
res(end+1) = ~isZero(C);

% no generator, but radius
C = capsule(zeros(4,1),zeros(4,1),1);
res(end+1) = ~isZero(C);

% not zero, but within tolerance
C = capsule(zeros(3,1),zeros(3,1),1);
tol = 2;
res(end+1) = isZero(C,tol);

% does not contain origin, but within tolerance
c = 0.5*[sqrt(2); sqrt(2)];
g = 0.5*[sqrt(2); -sqrt(2)];
r = 0.5;
tol = 2;
C = capsule(c,g,r);
res(end+1) = isZero(C,tol);

% combine tests
res = all(res);

%------------- END OF CODE --------------