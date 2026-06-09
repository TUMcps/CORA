function completed = example_polytope_halfspaces()
% example_polytope_halfspaces - show halfspaces of a polytope
%
% Syntax:
%    completed = example_polytope_halfspaces()
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

% Authors:       Tobias Ladner
% Written:       22-October-2015
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init polytope
A = [1 0 -1 0 1; 0 1 0 -1 1]';
b = [3; 2; 3; 2; 1];
P = polytope(A,b);

% visualize
figure; hold on;
plot(P,1:2,'DisplayName','Polytope'); enlargeAxis;
plotHalfspaces(P,1:2,'DisplayName','Halfspaces');
legend;

completed = true;

end

% ------------------------------ END OF CODE ------------------------------
