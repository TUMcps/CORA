function [R,tcomp] = observe_volMinA(obj,options)
% observe_volMinA - computes the guaranteed state estimation approach
% from [1]; the approach is extended here for nonlinear systems.
%
%
% Syntax:  
%    [R,tcomp] = observe_volMinA(obj,options)
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
%    [1] T. Alamo, J. M. Bravo, and E. F. Camacho. Guaranteed state 
%        estimation by zonotopes. Automatica, 41(6):1035–1043, 2005.
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
options.intersectionTechnique = 'alamo-volume';

% apply set-membership approach
tic
R = observe_stripBased(obj,options);
tcomp = toc;

end

%------------- END OF CODE --------------