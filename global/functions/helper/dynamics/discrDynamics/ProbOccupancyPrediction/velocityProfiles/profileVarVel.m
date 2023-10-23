function [velocity,input]=profileVarVel(pos,acc)
% profileVarVel - returns the velocity for a given position and maximum 
%    acceleration of the velocity profile of the corresponding path.
%
% Syntax:
%    [velocity,input] = profileVarVel(pos,acc)
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
% Written:       09-October-2009
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%speed limit:
velocity=20;

input=0;

% ------------------------------ END OF CODE ------------------------------
