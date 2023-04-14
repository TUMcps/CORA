function res = test_zonotope_mptPolytope
% test_zonotope_mptPolytope - unit test function of mptPolytope
%
% Syntax:  
%    res = test_zonotope_mptPolytope
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

% Author:       Matthias Althoff
% Written:      26-July-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% create zonotope
Z1 = zonotope([-4, -3, -2, -1; 1, 2, 3, 4]);

% obtain polytope
P = mptPolytope(Z1);

% extract MPT object
P = get(P,'P');

% obtain halfspace matrix
C = P.H(:,1:end-1);

% obtain distance vector
d = P.H(:,end);

% true results
true_C = [ 0.554700196225229   0.832050294337844; ...
           0.832050294337844   0.554700196225229; ...
           0.970142500145332   0.242535625036333; ...
          -0.554700196225229  -0.832050294337844; ...
          -0.832050294337844  -0.554700196225229; ...
          -0.970142500145332  -0.242535625036333];
      
true_d =  [2.773500981126146; 0; 0; 5.547001962252290; 5.547001962252290; 7.276068751089989];

% check result
res = compareMatrices([C,d],[true_C,true_d],1e-13);

%------------- END OF CODE --------------
