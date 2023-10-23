function res = test_equalDimCheck
% test_equalDimCheck - unit test function for dimensionality check
%
% Syntax:
%    res = test_equalDimCheck()
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

% Authors:       Mark Wetzlinger
% Written:       28-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

if ~CHECKS_ENABLED
    return
end

% set and set
Z = zonotope(zeros(2,1),eye(2));
I = interval(-ones(2,1),ones(2,1));
% check ok
equalDimCheck(Z,I);

% set and numeric vector
E = ellipsoid(eye(2),[1;-1]);
p = [2;1];
% check ok
equalDimCheck(E,p);

% square matrix and set
M = [2 3 1; -3 -2 0; 1 4 2];
I = interval([-3;5;1],[0;10;7]);
% check ok
equalDimCheck(M,I);

% non-square matrix and set
M = [2 1; -1 2; 0 3];
C = capsule([1;-1],[2;1],0.5);
% check ok
equalDimCheck(M,C);

% matrix set and set
intMat = intervalMatrix([-1 1; 0 1],[1 2; 1 0]);
pZ = polyZonotope([1;-1],[1 0 2; -1 1 -1],[0; 1],[2 1 0]);
% check ok
equalDimCheck(intMat,pZ);

% scalar and set
s = 5;
Z = zonotope([1;-1;4],[1 9 3 -4; 2 0 1 -2; 8 6 -1 8]);
% check ok
equalDimCheck(s,Z);


% dimension mismatches

% set and set
Z = zonotope(zeros(3,1),eye(3));
I = interval(-ones(2,1),ones(2,1));
try 
    equalDimCheck(Z,I);
    res = false;
end

% set and numeric vector
E = ellipsoid(eye(2),[1;-1]);
p = [2;1;-1];
try 
    equalDimCheck(E,p);
    res = false;
end

% square matrix and set
M = [2 3 1; -3 -2 0; 1 4 2];
I = interval([-3;5],[0;10]);
try
    equalDimCheck(M,I);
    res = false;
end

% non-square matrix and set
M = [2 1; -1 2; 0 3];
C = capsule([1;-1;1],[2;1;-1],0.5);
try
    equalDimCheck(M,C);
    res = false;
end

% matrix set and set
intMat = intervalMatrix([-1 1 2; 0 1 1],[1 2 0; 1 0 3]);
pZ = polyZonotope([1;-1],[1 0 2; -1 1 -1],[0; 1],[2 1 0]);
try
    equalDimCheck(intMat,pZ);
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
