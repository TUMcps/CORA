function coordinateMat = gridPoints(obj,segments)
% gridPoints - Computes grid points of an interval; the points are
% generated in a way, such that a continuous space is uniformly partitioned.
%
% Syntax:  
%    coordinateMat = gridPoints(obj,segments)
%
% Inputs:
%    obj - interval object
%    segments - number of segments for each dimension (scalar of column vector)
%
% Outputs:
%    coordinateMat - matrix where columns are the grid points
%
% Example: 
%    I = interval([1 2; -1 1]);
%    coordinates = gridPoints(I,10);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      22-July-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%obtain segment length
segLengthVec = (obj.sup-obj.inf)./segments;

%obtain first grid point
startingPoint = obj.inf(:,1) + 0.5*segLengthVec;
coordinateMat = startingPoint;

%obtain segment combinations
dim = length(obj.inf(:,1));
if isscalar(segments)
    comb = full_fact_mod(ones(dim,1)*segments);
else
    comb = full_fact_mod(segments);
end

%add further grid points
for i = 2:length(comb(:,1))
    coordinateMat(:,i) = startingPoint+(comb(i,:)'-1).*segLengthVec;
end

%------------- END OF CODE --------------