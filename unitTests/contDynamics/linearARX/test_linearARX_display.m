function res = test_linearARX_display
% test_linearARX_display - unit test for computing a GO model for a
%   linearARX object
%
% Syntax:
%    res = test_linearARX_display
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Laura Luetzow
% Written:       05-July-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% initialize system and display
A_bar = {[-0.4 0.6; 0.6 -0.4];[0.1 0; 0.2 -0.5]};
B_bar = {[0; 0];[0.3; -0.7];[0.1; 0]};
dt = 0.1;
sys = linearARX('exampleSys',A_bar, B_bar,dt)

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
