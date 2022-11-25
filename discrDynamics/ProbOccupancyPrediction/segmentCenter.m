function [c] = segmentCenter(segments, segVector, devVector, widthOfLane, nrOfDevSegments)
% segmentCenter - returns the centers of reachable set segments where the
% segment number vector(along the path) and the deviation number vector
% (deviation from the path) is given
%
% replaces segCenter.m for using without the CORA road class

%obtain x,y and angle values
x = segments.x;
y = segments.y;
angle = segments.angle;

%get segment centers along the path
cSeg = [];
for i = 1:length(segVector)
    %obtain segment index
    ind = segVector(i);
    %compute segment center
    cSeg(end+1,:) = [x(ind), y(ind)] + 0.5 * [x(ind+1) - x(ind), y(ind+1) - y(ind)];
end

%check if additionally deviation should be considered
if isempty(devVector)
    c = cSeg;
else
    %compute centers due to deviation in addition
    cDev = [];
    for i = 1:length(segVector)
        %obtain segment index
        ind = segVector(i);
        
        %compute mean angle
        angle = angle(ind) + 0.5 * (angle(ind+1) - angle(ind));
        
        %calculate translation vector tangential to the path
        transLat(1,1) = cos(angle - 0.5 * pi) * widthOfLane / nrOfDevSegments;
        transLat(1,2) = sin(angle - 0.5 * pi) * widthOfLane / nrOfDevSegments; 
        
        %get center
        cTmp = cSeg(i,:);
        %deviation loop
        for iDev = 1:nrOfDevSegments
            %obtain deviation index
            devInd = devVector(iDev);
            %compute deviation center
            cDev(end+1,:) = cTmp + (devInd - 0.5 - 0.5 * nrOfDevSegments) * transLat;
        end
    end
    c = cDev;
end