function [velocity,input]=profile4(pos,acc)
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
% Written:      02-July-2008 
% Last update:  09-October-2009
% Last revision:---

%------------- BEGIN CODE --------------

%compute absolute value of acceleration
acc=abs(acc);

%set speed limit
speedLimit=16; %[m/s]
radius=15.5; %[m]
curveSpeed=sqrt(acc*radius); %[m/s]

deltaT=(speedLimit-curveSpeed)/acc;
deltaX=0.5*acc*deltaT^2+curveSpeed*deltaT;
pos1=60-deltaX;
pos2=60;
pos3=85;
pos4=85+deltaX;
pos5=inf;

%determine position segments
posSegments=[pos1,pos2,pos3,pos4,pos5];

%find position segment
ind = find(pos<=posSegments, 1, 'first');

%select position segment
switch ind
  case 1
    velocity=speedLimit;
    input=0;
  case 2
    velocity=sqrt(speedLimit^2-2*(pos-pos1)*acc);
    input=-1;
  case 3
    velocity=sqrt(acc*radius);
    input=0;
  case 4
    velocity=sqrt(curveSpeed^2+2*(pos-pos3)*acc);
    input=1;
  case 5
    velocity=speedLimit;
    input=0;
  otherwise
    disp('Velocity profile error');
end

%------------- END OF CODE --------------