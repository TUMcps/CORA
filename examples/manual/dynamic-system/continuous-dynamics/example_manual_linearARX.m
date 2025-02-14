function example_manual_linearARX()
% example_manual_linearARX - example from the manual demonstrating the 
% linearARX constructor as defined in the manual
%
% Syntax:
%   example_manual_linearARX()
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

% Authors:       Laura Luetzow
% Written:       11-October-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% system matrices
A_bar = {[0.4 0.6; 0.6 0.4];
    [0.1 0; 0.2 0.1]};
B_bar = {[0; 0];[0.3; -0.7];
    [0.1; 0]};

% sampling time
dt = 0.1;

% ARX system
sys = linearARX(A_bar, B_bar, dt);

% ------------------------------ END OF CODE ------------------------------
