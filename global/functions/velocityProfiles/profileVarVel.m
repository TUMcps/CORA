function [velocity,input]=profileVarVel(pos,acc)
% profile1..n - returns the velocity for a given position and maximum 
% accelerationof the velocity profile of the corresponding path.
%
% Syntax:  
%    [velocity]=profile1(pos)
%
% Inputs:
%    pos - position on the path
%
% Outputs:
%    velocity - velocity of the velocity profile
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      09-October-2009
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%speed limit:
velocity=20;

input=0;


%------------- END OF CODE --------------