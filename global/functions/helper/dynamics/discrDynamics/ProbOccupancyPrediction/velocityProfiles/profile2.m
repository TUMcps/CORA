function [velocity,input] = profile2(pos,acc)
% profile2 - returns the velocity for a given position and maximum 
%    acceleration of the velocity profile of the corresponding path.
%
% Syntax:
%    [velocity,input] = profile2(pos,acc)
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

%path 1 has a speed limit of 15
velocity=16;

input=0;

% ------------------------------ END OF CODE ------------------------------
