function [velocity,input]=profilePlot(pos,acc)
% profilePlot - returns the velocity for a given position and maximum 
% acceleration of the velocity profile of the corresponding path.
%
% Syntax:
%    [velocity,input]=profilePlot(pos,acc)
%
% Inputs:
%    pos - position on the path
%    acc - ???
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

%set speed limit
speedLimit=16.6; %[m/s]
radius=15.5; %[m]
curveSpeed=sqrt(acc*radius); %[m/s]

deltaT=(speedLimit-curveSpeed)/acc;
deltaX=0.5*acc*deltaT^2+curveSpeed*deltaT;
pos1=65-deltaX;
pos2=65;
pos3=90;
pos4=90+deltaX;
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

% ------------------------------ END OF CODE ------------------------------
