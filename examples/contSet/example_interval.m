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

I1 = interval([0; -1], [3; 1]); % create interval I1
I2 = interval([-1; -1.5], [1; -0.5]); % create interval I2

r = rad(I1) % obtain and display radius of I1
I3 = I1 & I2; % computes the intersection of I1 and I2
c = center(I3) % returns and displays the center of I3

figure; hold on;
plot(I3,[1 2],'FaceColor',colorblind('b')); % plot I3 
plot(I1); % plot I1
plot(I2); % plot I2

%example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
