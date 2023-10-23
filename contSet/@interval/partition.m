function dzNew = partition(I, splits)
% partition - partitions a multidimensional interval into subintervals
%
% Syntax:
%    I = partition(I, splits)
%
% Inputs:
%    I - interval object
%    splits - number of splits in each dimension
%
% Outputs:
%    dzNew - cell array of intervals
%
% Example:
%    I = interval([2;3;4],[5;6;7]);
%    Isplit = partition(I,2);
%    figure;
%    for i=1:length(Isplit)
%        subplot(1,2,1); hold on;
%        plot(Isplit{i},[1,2]);
%        subplot(1,2,2); hold on;
%        plot(Isplit{i},[2,3]);
%    end
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       19-September-2012
% Last update:   04-May-2020 (MW, migrated "splitIntervals" to @interval)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input arguments
inputArgsCheck({{I,'att','interval'};
                {splits,'att','numeric',{'nonnan','scalar','positive'}}});

%increment
incr = (supremum(I) - infimum(I)) / splits;
inf = infimum(I);

%possible indices
ind = combinator(splits,dim(I),'p','r');

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


% ------------------------------ END OF CODE ------------------------------
