function dzNew = partition(Int, splits)
% partition - partitions a multidimensional interval into subintervals
%
% Syntax:  
%    Int = partition(Int, splits)
%
% Inputs:
%    Int - interval vector
%    splits - number of splits in each dimension
%
% Outputs:
%    dzNew - cell array of intervals
%
% Example:
%    I = interval([2;3;4],[5;6;7]);
%    Isplit = partition(I,2);
%    figure; hold on;
%    for i=1:length(Isplit)
%       subplot(1,2,1); hold on;
%       plot(Isplit{i},[1,2]);
%       subplot(1,2,2); hold on;
%       plot(Isplit{i},[2,3]);
%    end
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      19-September-2012
% Last update:  04-May-2020 (MW, migrated "splitIntervals" to @interval)
% Last revision:---

%------------- BEGIN CODE --------------

%increment
incr = (supremum(Int) - infimum(Int)) / splits;
inf = infimum(Int);

%possible indices
ind = combinator(splits,dim(Int),'p','r');

nrOfComb = length(ind(:,1));
dzNew = cell(nrOfComb,1);
%loop
for i = 1:nrOfComb
    %get current indices
    currInd = ind(i,:);
    %bounds
    lowerBound = inf + (currInd'-1).*incr;
    upperBound = inf + (currInd').*incr;
    %new interval
    dzNew{i} = interval(lowerBound, upperBound);
end


%------------- END OF CODE --------------