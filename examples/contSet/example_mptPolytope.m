function completed = example_mptPolytope()
% example_mptPolytope - example instantiation of mptPolytope objects
%
% Syntax:  
%    completed = example_mptPolytope()
%
% Inputs:
%    -
%
% Outputs:
%    completed - boolean
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% Author:        Matthias Althoff
% Written:       21-April-2018
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

Z1 = zonotope([1 1 1; 1 -1 1]); % create zonotope Z1
Z2 = zonotope([-1 1 0; 1 0 1]); % create zonotope Z2

P1 = polytope(Z1); % convert zonotope Z1 to halfspace representation
P2 = polytope(Z2); % convert zonotope Z2 to halfspace representation

P3 = P1 + P2 % perform Minkowski addition and display result
P4 = P1 & P2; % compute intersection of P1 and P2

V = vertices(P4) % obtain and display vertices of P4

figure; hold on
plot(P1); % plot P1
plot(P2); % plot P2
plot(P3,[1 2],'g'); % plot P3
plot(P4,[1 2],'FaceColor',[.6 .6 .6],'Filled',true,'EdgeColor','none'); % plot P4 

%example completed
completed = 1;

%------------- END OF CODE --------------