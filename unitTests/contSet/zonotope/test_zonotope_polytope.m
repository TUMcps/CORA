function res = test_zonotope_polytope
% test_zonotope_polytope - unit test function of polytope
%
% Syntax:  
%    res = test_zonotope_polytope
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

% Author:       Matthias Althoff
% Written:      26-July-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

% Analytical Test ---------------------------------------------------------

% create zonotope
Z1 = zonotope([-4, -3, -2, -1; 1, 2, 3, 4]);

% obtain polytope
P = polytope(Z1);

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
for i = size(C,1)
    [temp,ind] = ismembertol(C(i,:),true_C,1e-13,'ByRows',true);
    if ~temp || abs(d(i)-true_d(ind)) > 1e-13 || size(C,1) ~= size(true_C,1)
        error('Analytical test failed!'); 
    end
end

if res
    disp('test_zonotope_polytope successful');
end

%------------- END OF CODE --------------
