function dzNew = splitIntervals(dz, splits)
% splitIntervals - splits a multidimensional interval into subintervals
%
% Syntax:  
%    dz = splitIntervals(dz, splits)
%
% Inputs:
%    dz - interval vector
%    splits - number of splits in each dimension
%
% Outputs:
%    dzNew - cell array of intervals
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      19-September-2012
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%get dimension
dim = length(dz);

%increment
inc = (supremum(dz) - infimum(dz))/splits;
inf = infimum(dz);

%possible indices
ind = combinator(splits,dim,'p','r');

%loop
for i = 1 : length(ind(:,1))
    %get current indices
    currInd = ind(i,:);
    %bounds
    lowerBound = inf + (currInd'-1).*inc;
    upperBound = inf + (currInd').*inc;
    %new interval
    dzNew{i} = interval(lowerBound, upperBound);
end



%------------- END OF CODE --------------