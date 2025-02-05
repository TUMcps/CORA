function plotPositionDistribution(segments, segmentProb, devProb, widthOfLane)
% plotPositionDistribution - replaces plot.m in @road without using the CORA road class
%
% Syntax:
%    plotPositionDistribution(segments, segmentProb, devProb, widthOfLane)
%
% Inputs:
%    ???
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       ???
% Written:       ???
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%get number of deviation segments
nrOfDev = length(devProb);

%find segments that should be plotted => all segments with a probability > 0
ind = find(segmentProb);

for iInd = 1:length(ind)
    %get segment number
    iSeg = ind(iInd);
    if segmentProb(iSeg) > 0.0001
        %iDev: deviation segment
        for iDev = 1:nrOfDev
            %obtain segment polytope
            [P] = getPolytopeForSegment(segments, iSeg, iDev, nrOfDev, widthOfLane);
            
            %plot segment
            %for more options on plotting see the origin file
            V = vertices(extreme(P)');
            plot(V, 'grayTones', segmentProb(iSeg) * devProb(iDev));

            hold on
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
