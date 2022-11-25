function [R,tcomp] = observe_volMinB(obj,options)
% observe_volMinA - computes the guaranted state estimation approach
% from [1]; the approach is extended here for nonlinear systems.
%
%
% Syntax:  
%    [R,tcomp] = observe_volMinB(obj,options)
%
% Inputs:
%    obj - discrete-time nonlinear system object
%    options - options for the guaranteed state estimation
%
% Outputs:
%    R - observed set of points in time
%    tcomp - computation time
%
% Reference:
%    [1] J. M. Bravo, T. Alamo, and E. F. Camacho. Bounded error
%        identification of systems with time-varying parameters. IEEE
%        Transactions on Automatic Control, 51(7):1144–1150, 2006.
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Matthias Althoff
% Written:       27-Mar-2021
% Last update:   ---
% Last revision: ---


%------------- BEGIN CODE --------------

% set intersection procedure
options.intersectionType = 1;
options.intersectionTechnique = 'bravo'; % unclear whether method is an intersection method

% apply set-membership approach
tic
R = observe_stripBased(obj,options);
tcomp = toc;

end

%------------- END OF CODE --------------