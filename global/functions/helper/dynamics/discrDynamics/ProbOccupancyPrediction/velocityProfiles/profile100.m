function [velocity,input] = profile100(pos,acc)
% profile100 - returns the velocity for a given position and maximum 
%    acceleration of the velocity profile of the corresponding path.
%
% Syntax:
%    [velocity,input] = profile100(pos,acc)
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
velocity=100/3.6;
%velocity=7;

input=0;

% ------------------------------ END OF CODE ------------------------------
