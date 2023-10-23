function acceleration = accelerationConstraints(velocity,acceleration,p)
% accelerationConstraints - adjusts the acceleration based on acceleration
% constraints
%
% Syntax:  
%    acceleration = accelerationConstraints(velocity,acceleration,p)
%
% Inputs:
%    acceleration - acceleration in driving direction
%    velocity - velocity in driving direction
%    p - longitudinal parameter structure
%
% Outputs:
%    acceleration - acceleration in driving direction
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      15-December-2017
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%positive acceleration limit
if velocity>p.v_switch
    posLimit = p.a_max*p.v_switch/velocity;
else
    posLimit = p.a_max;
end

%acceleration limit reached?
if (velocity<=p.v_min && acceleration<=0) || (velocity>=p.v_max && acceleration>=0)
    acceleration = 0;
elseif acceleration<=-p.a_max
    acceleration = -p.a_max;
elseif acceleration>=posLimit
    acceleration = posLimit; 
end

%------------- END OF CODE --------------
