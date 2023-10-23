function res = test_interval_capsule
% test_interval_capsule - unit test function of conversion to capsules
%
% Syntax:
%    res = test_interval_capsule
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
% Written:       28-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init interval
I = interval([-2;-1],[4;7]);

% convert to capsule
C = capsule(I);
% convert back to interval
II = interval(C);

% check containment of vertices of interval
V = vertices(I);
res = all(contains(C,V));
% check containment of interval
res(end+1,1) = contains(II,I);


% degenerate interval
I = interval([-2;0;5],[-1;0;7]);

% convert to capsule
C = capsule(I);

% check containment of vertices of interval
V = vertices(I);
res(end+1,1) = all(contains(C,V));


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
