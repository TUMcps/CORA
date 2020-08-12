function [velocity,input]=profileNoLimit(pos,acc)
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

%path 1 has a speed limit of 1e3
velocity=1e3;

input=0;

%------------- END OF CODE --------------