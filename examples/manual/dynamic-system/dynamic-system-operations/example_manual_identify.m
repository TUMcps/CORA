function example_manual_identify()
% example_manual_identify - example from the manual demonstrating the 
% identify operation as defined in the manual
%
% Syntax:
%   example_manual_identify()
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

% Authors:       Niklas Kochdumper
% Written:       15-December-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% original system
A = [-0.7 -2; 2 -0.7];
sysOrig = linearSys(A);

% simulation of the original system
simOpts.x0 = [10; 5];
simOpts.tFinal = 5;
[t,x] = simulate(sysOrig,simOpts);

% system identification
sys = linearSys.identify(x,t)

% ------------------------------ END OF CODE ------------------------------
