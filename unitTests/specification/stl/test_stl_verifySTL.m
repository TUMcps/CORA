function res = test_stl_verifySTL
% test_stl_verifySTL - unit test function of verifySTL
%
% Syntax:
%    res = test_stl_verifySTL
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
% See also: stl

% Authors:       Benedikt Seidl
% Written:       15-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% specification
x = stl('x',2);
eq = until(x(2) > 0.4,x(1) < 0,interval(0,2));

% dynamic system
sys = linearSys([0 -1; 1 0],[0; 0]);

% reachability parameters
params.R0 = zonotope(interval([0.5;0.5],[1;1]));

% algorithm parameters
options.verifyAlg = 'stl:seidl';
options.taylorTerms = 10;
options.zonotopeOrder = 10;

% verify
res = verify(sys,params,options,eq);

% ------------------------------ END OF CODE ------------------------------
