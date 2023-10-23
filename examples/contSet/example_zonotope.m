function completed = example_zonotope()
% example_zonotope - example instantiation of zonotope objects
%
% Syntax:
%    completed = example_zonotope()
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
% Written:       21-April-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

Z1 = zonotope([1 1 1; 1 -1 1]); % create zonotope Z1
Z2 = zonotope([-1 1 0; 1 0 1]); % create zonotope Z2
A = [0.5 1; 1 0.5]; % numerical matrix A

Z3 = Z1 + Z2; % Minkowski addition
Z4 = A*Z3; % linear map

figure; hold on
plot(Z1,[1 2],'b'); % plot Z1 in blue
plot(Z2,[1 2],'g'); % plot Z2 in green
plot(Z3,[1 2],'r'); % plot Z3 in red
plot(Z4,[1 2],'k'); % plot Z4 in black

P = polytope(Z4) % convert to and display halfspace representation
I = interval(Z4) % convert to and display interval

figure; hold on
plot(Z4); % plot Z4
plot(I,[1 2],'g'); % plot interval in green

%example completed
completed = 1;

% ------------------------------ END OF CODE ------------------------------
