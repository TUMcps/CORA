function [iS, iProp] = exactIntersectingSegments(Obj,contSet)
% exactIntersectingSegments - finds the exact segments of the partition 
% that intersect a set P, and the proportion of P that is in each segment.
%
% Syntax:  
%    [iS, iProp] = exactIntersectingSegments(Obj,contSet)
%
% Inputs:
%    Obj - partition object
%    contSet - continuous set
%
% Outputs:
%    iS - indizes of the partitions that intersect the set P (index starts
%         at 0)
%    iProp - proportion of the volume of P in the corresponding segment in
%            iS (including the zero element = environment outside
%            partition)
%
% Example: 
%    -
%
% Other m-files required: tbd
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Matthias Althoff, Aaron Pereira, Niklas Kochdumper
% Written:      09-October-2006 (MA)
% Last update:  02-August-2017 (AP)
%               13-November-2017 (MA)
%               26-March-2018 (NK)
% Last revision:---

%------------- BEGIN CODE --------------

%% change the contSet to a polytope

if isa(contSet,'mptPolytope')
    P = contSet;
else
    P = polytope(contSet);
end


% initialize
iS = [];
iProp = [];

% compute volume of P
volTotal = volume(P);

% if P is not a singleton
if volTotal > 0 
    % obtain cells that intersect with P
    cellIndices = intersectingCells(Obj,P);
    polytopes_to_check = cellPolytopes(Obj,cellIndices);
    
    % reserve memory for speed
    nrOfCells = prod(Obj.nrOfSegments);
    iProp = zeros(nrOfCells,1);
    
    % compute intersection and volume for each cell
    for i = 1:length(polytopes_to_check)
        sP = polytopes_to_check{i};
        iP = P&sP; %iP:intersected polytope
        volPartial = volume(iP);
        if volPartial > 0 
            iS = [iS;cellIndices(i)];
            iProp(iS(end),1) = volPartial/volTotal;
        end
    end
end

% epsilon value for numeric imprecision
epsilon = 1e-6;

% some of the volume might be outside the partition:
proportionOutside = 1-sum(iProp);
% if the proportions sum to more than one, something is wrong!
if proportionOutside < -epsilon
    disp('(transition) probability error')
    disp(['P(outside)=',num2str(proportionOutside)])
end

% account for "zero" segment
iS = [0; iS];
iProp = [proportionOutside; iProp];

%------------- END OF CODE --------------