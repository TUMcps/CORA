function completed = example_interval()
% example_interval - example instantiation of interval objects
%
% Syntax:  
%    completed = example_interval()
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

I1 = interval([0; -1], [3; 1]); % create interval I1
I2 = interval([-1; -1.5], [1; -0.5]); % create interval I2
Z1 = zonotope([1 1 1; 1 -1 1]); % create zonotope Z1

r = rad(I1) % obtain and display radius of I1
is_intersecting = isIntersecting(I1, Z1) % Z1 intersecting I1?
I3 = I1 & I2; % computes the intersection of I1 and I2
c = center(I3) % returns and displays the center of I3

figure; hold on
plot(I1); % plot I1
plot(I2); % plot I2
plot(Z1,[1 2],'g'); % plot Z1
plot(I3,[1 2],'FaceColor',[.6 .6 .6],'Filled',true,'EdgeColor','none'); % plot I3 

%example completed
completed = 1;

%------------- END OF CODE --------------