function [c]=segCenter(varargin)
% segCenter - returns the centers of reachable set segments where the
% segment number vector(along the path) and the deviation number vector
% (deviation from the path) is given
%
% Syntax:  
%    [c]=segCenter(obj,segVector,devVector)
%
% Inputs:
%    obj - road object
%    segVector - vector of indices along the path
%    devVector - vector of indices orthogonal to the path 
%
% Outputs:
%    c - segment centers
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

%get inputs
if nargin==2
    obj=varargin{1};
    segVector=varargin{2};
    devVector=[];
elseif nargin==3
    obj=varargin{1};
    segVector=varargin{2};
    devVector=varargin{3};    
end

%obtain x,y and angle values
x=obj.segments.x;
y=obj.segments.y;
angle=obj.segments.angle;

%get segment centers along the path
cSeg=[];
for i=1:length(segVector)
    %obtain segment index
    ind=segVector(i);
    %compute segment center
    cSeg(end+1,:)=[x(ind),y(ind)]+0.5*[x(ind+1)-x(ind),y(ind+1)-y(ind)];
end

%check if additionally deviation should be considered
if isempty(devVector)
    c=cSeg;
else
    %compute centers due to deviation in addition
    cDev=[];
    for i=1:length(segVector)
        %obtain segment index
        ind=segVector(i);
        
        %compute mean angle
        angle=angle(ind)+0.5*(angle(ind+1)-angle(ind));
        %calculate translation vector tangential to the path
        transLat(1,1)=cos(angle-0.5*pi)*obj.width/obj.nrOfDevSegments;
        transLat(1,2)=sin(angle-0.5*pi)*obj.width/obj.nrOfDevSegments; 
        
        %get center
        cTmp=cSeg(i,:);
        %deviation loop
        for iDev=1:length(devVector)
            %obtain deviation index
            devInd=devVector(iDev);
            %compute deviation center
            cDev(end+1,:)=cTmp+(devInd-0.5-0.5*obj.nrOfDevSegments)*transLat;
        end
    end
    c=cDev;
end
    

%------------- END OF CODE --------------