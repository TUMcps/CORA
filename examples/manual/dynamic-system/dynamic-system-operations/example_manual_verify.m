function example_manual_verify()
% example_manual_verify - example from the manual demonstrating the 
% verify operation as defined in the manual
%
% Syntax:
%   example_manual_verify()
%
% Inputs:
%    -
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Authors:       Benedikt Seidl, Tobias Ladner
% Written:       19-October-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% specification
x = stl('x',2);
spec = until(x(2) > 0.4,x(1) < 0,interval(0,2));

% dynamic system
sys = linearSys([0 -1; 1 0],[0; 0]);

% reachability parameters
params.R0 = zonotope(interval([0.5;0.5],[1;1]));

% algorithm settings
options.verifyAlg = 'stl:seidl';
options.taylorTerms = 10;
options.zonotopeOrder = 10;

% verify
verify(sys,params,options,spec)

% ------------------------------ END OF CODE ------------------------------
