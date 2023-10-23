function completed = example_taylm()
% example_taylm - example instantiation of taylm objects
%
% Syntax:
%    completed = example_taylm()
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       29-June-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

a1 = interval(-1, 2);   % generate a scalar interval [-1,2]
a2 = interval(2, 3);    % generate a scalar interval [2,3]
a3 = interval(-6, -4);  % generate a scalar interval [-1,2]
a4 = interval(4, 6);    % generate a scalar interval [2,3]

b1 = taylm(a1, 6);      % Taylor model with maximum order of 6 and name a1
b2 = taylm(a2, 6);      % Taylor model with maximum order of 6 and name a2
b3 = taylm(a3, 6);      % Taylor model with maximum order of 6 and name a3
b4 = taylm(a4, 6);      % Taylor model with maximum order of 6 and name a4

B1 = [b1; b2]   % generate a row of Taylor models
B2 = [b3; b4]   % generate a row of Taylor models

B1 + B2     % addition
B1' * B2    % matrix multiplication
B1 .* B2    % pointwise multiplication
B1/2        % division by scalar
B1./B2      % pointwise division
B1.^3       % power function
sin(B1)     % sine function

sin(B1(1,1)) + B1(2,1).^2 - B1' * B2 % combination of functions


completed = true;

% ------------------------------ END OF CODE ------------------------------
