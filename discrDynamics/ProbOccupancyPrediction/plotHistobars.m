function plotHistobars(obj,posProb)
% plotHistobars - plots the a probability distribution, similar to
% plotHisto
%
% Syntax:  
%    plotHistobars(obj,posProb)
%
% Inputs:
%    obj - partition object
%    posProb - position probability vector
%
% Outputs:
%    ---
%
% Example: 
%    Text for example...
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Matthias Althoff
% Written:      05-March-2009
% Last update:  25-July-2016
% Last revision:---

%------------- BEGIN CODE --------------

%get position interval length
posInt=obj.nrOfSegments(1);

%generate position vector
posVec=linspace(obj.intervals(1,1),obj.intervals(1,2),posInt+1);


%figure;

%init x,y
x(1)=posVec(1);
y(1)=posProb(1);
x(2)=posVec(2);
y(2)=y(1);

%loop
for i=1:obj.nrOfSegments(1)
   %instantiate interval
   IH=interval([posVec(i);0],[posVec(i+1);posProb(i)]);
   %plot interval hull
   plot(IH,[1 2],'grayFrame');
end


%------------- END OF CODE --------------