function completed = example_zonoBundle()
% example_zonoBundle - example instantiation of zonoBundle objects
%
% Syntax:
%    completed = example_zonoBundle()
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

Z{1} = zonotope([1 1 1; 1 -1 1]); % create zonotope Z1;
Z{2} = zonotope([-1 1 0; 1 0 1]); % create zonotope Z2;
Zb = zonoBundle(Z); % instantiate zonotope bundle from Z1, Z2
vol = volume(Zb) % compute and display volume of zonotope bundle

figure; hold on
plot(Zb,[1 2],'FaceColor',colorblind('gray')); % plot Zb in gray
plot(Z{1}); % plot Z1 
plot(Z{2}); % plot Z2 

%example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
