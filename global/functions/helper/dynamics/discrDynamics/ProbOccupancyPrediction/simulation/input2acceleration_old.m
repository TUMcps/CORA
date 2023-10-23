function acceleration = input2acceleration_old(input,velocity,type)
% input2acceleration_old - transforms the input (-1,1) to a physical 
% acceleration. For decelartion with input (-1,0), the deceleration is
% modeled linear with a maximum deceleration of -10 m/s^2. For the
% acceleration, some linearization of the acceleration function has to be
% performed.
%
% Syntax:
%    acceleration = input2acceleration_old(input,velocity,type)
%
% Inputs:
%    input - input of the car ranging from (-1,1)
%    velocity - velocity of the car
%
% Outputs:
%    acceleration - acceleration value of the car
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       30-June-2008
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%set constants for acceleration model
if strcmp(type,'car')
    c4=10;
    c5=60;
elseif strcmp(type,'bicycle')
    c4=5;
    c5=7;
end

%if the car brakes
if input<0
    acceleration=10*input; %max deceleration: 10 m/s^2
%if the car accelerates
else
    acceleration=c4*(1-((velocity/c5))^2)*input; %speed
end

% ------------------------------ END OF CODE ------------------------------
