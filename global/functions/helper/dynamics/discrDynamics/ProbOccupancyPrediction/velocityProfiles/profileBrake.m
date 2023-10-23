function [velocity,input] = profileBrake(pos,acc)
% profileBrake - returns the velocity for a given position and maximum 
%    acceleration of the velocity profile of the corresponding path.
%
% Syntax:
%    [velocity,input] = profileBrake(pos,acc)
%
% Inputs:
%    pos - position on the path
%    acc - acceleration
%
% Outputs:
%    velocity - velocity of the velocity profile
%    input - ???
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       02-July-2008 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

initSpeed=12;

%determine velocity
velocity=sqrt(initSpeed^2-2*(pos)*acc);

% ------------------------------ END OF CODE ------------------------------
