function res = test_mptPolytope_generateRandom
% test_mptPolytope_generateRandom - unit test function of generateRandom
%
% Syntax:  
%    res = test_mptPolytope_generateRandom
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
% Written:      19-May-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% empty call
P = mptPolytope.generateRandom();

% values for tests
n = 3;
c = [2;1;-1];

% only dimension
P = mptPolytope.generateRandom('Dimension',n);
res = dim(P) == n;

% only center (no check)
P = mptPolytope.generateRandom('Center',c);
res(end+1,1) = dim(P) == length(c);

% dimension and center
P = mptPolytope.generateRandom('Dimension',n,'Center',c);
res(end+1,1) = dim(P) == n;


% unify results
res = all(res);

%------------- END OF CODE --------------