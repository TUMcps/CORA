function [velocity,input]=profileBrake(pos,acc)
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

% Author: Matthias Althoff
% Written: 02-July-2008 
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

initSpeed=12;

%determine velocity
velocity=sqrt(initSpeed^2-2*(pos)*acc);


%------------- END OF CODE --------------