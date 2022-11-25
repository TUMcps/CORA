function plotVehSeg(obj,iSeg,iDevSeg,Zvehicle)
% plotVehSeg - plots the region a vehicle possibly occupies when given the
% segment index, the deviation index and the angle of the vehicle.
%
% Syntax:  
%    plotVehSeg(obj,iSeg,iDevSeg,Zvehicle)
%
% Inputs:
%    obj - road object
%    iSeg - path segment
%    iDevSeg - deviation segment
%    Zvehicle - unoriented and untranslated region of the vehicle specified
%    as zonotope
%
% Outputs:
%    ---
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 28-January-2009
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------


%get center and angle of the vehicle
[c,angle]=segData(obj,iSeg,iDevSeg);

%generate rotation matrix
Rot=[cos(angle) -sin(angle);...
     sin(angle) cos(angle)];

%rotate and translate vehicle region
Zmapped=Rot*Zvehicle+c;

%plot region
plot(Zmapped);


%------------- END OF CODE --------------