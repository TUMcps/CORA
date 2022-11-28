function res = test_halfspace_generateRandom
% test_halfspace_generateRandom - unit test function of generateRandom
%
% Syntax:  
%    res = test_halfspace_generateRandom
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
% Written:      27-Sep-2019
% Last update:  19-May-2022 (name-value pair syntax)
% Last revision:---

%------------- BEGIN CODE --------------

% empty call
h = halfspace.generateRandom();

% values for tests
n = 3;
c = [2;1;-1];
d = 3.5;

% only dimension
h = halfspace.generateRandom('Dimension',n);
res = dim(h) == n;

% only normal vector
h = halfspace.generateRandom('NormalVector',c);
res(end+1,1) = all(abs(h.c - c) < eps);

% only offset
h = halfspace.generateRandom('Offset',d);
res(end+1,1) = abs(h.d - d) < eps;

% dimension and normal vector
h = halfspace.generateRandom('Dimension',n,'NormalVector',c);
res(end+1,1) = dim(h) == n && all(abs(h.c - c) < eps);

% dimension and offset
h = halfspace.generateRandom('Dimension',n,'Offset',d);
res(end+1,1) = dim(h) == n && abs(h.d - d) < eps;

% normal vector and offset
h = halfspace.generateRandom('NormalVector',c,'Offset',d);
res(end+1,1) = all(abs(h.c - c) < eps) && abs(h.d - d) < eps;

% dimension, normal vector, and offset
h = halfspace.generateRandom('Dimension',n,'NormalVector',c,'Offset',d);
res(end+1,1) = dim(h) == n && all(abs(h.c - c) < eps) && abs(h.d - d) < eps;


% unify results
res = all(res);

%------------- END OF CODE --------------