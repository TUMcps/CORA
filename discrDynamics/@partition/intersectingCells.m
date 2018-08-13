function [cells,error] = intersectingCells(obj,contSet,varargin)
% intersectingCells - returns the cells possibly intersecting with a
% continuous set, overapproximatively, by overapproximating the convex set 
% as a multidimensional interval.
%
% Syntax:  
%   [cells,error] = intersectingCells(obj,contSet,varargin)
%
% Inputs:
%    obj - partition object
%    contSet - a continuous set
%    varargin{1} - 'subscripts' gives answer as subscripts of the 
%                   partition cells, 'indices' gives the answer as 
%                   indices. Default is indices.
%
% Outputs:
%    cells - either the indices or the subscripts of the possibly intersecting cells
%    error - error flag
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff, Aaron Pereira
% Written:      14-September-2006
% Last update:  16-August-2007 
%               29-October-2007
%               17-September-2015
%               02-August-2017 (AP)
%               08-August-2018 (MA)
% Last revision:---


%------------- BEGIN CODE --------------

if nargin > 2
    switch varargin{1}
        case 'subscripts'
            giveAnswerAsIndices = 0;
        case 'indices'
            giveAnswerAsIndices = 1;
        otherwise
            disp('invalid choice of output, defaulting to indices')
            giveAnswerAsIndices = 1;
    end
else
    % default case is indices
    giveAnswerAsIndices = 1;
end


if isa(contSet,'zonotope')
    I=interval(contSet);
elseif isa(contSet,'polytope')
    V=vertices(extreme(contSet)');
    I=interval(V);    
elseif isa(contSet,'mptPolytope')
    I=interval(contSet);  
else
    I=contSet;
end


error=0; % set error to 0
if isa(I, 'interval')
    leftLimit = infimum(I);
    rightLimit = supremum(I);

else
    leftLimit = I;
    rightLimit = I;
end

if (size(I,2) ~= 1)
    if (size(I,1) ~= 1)
        disp('either define an interval or a point')
        return
    else
        leftLimit = leftLimit';
        rightLimit = rightLimit';
    end
end

bounds = [leftLimit,rightLimit];

if size(bounds,1) ~= length(obj.nrOfSegments)
    disp('state space and subset are not of same dimension')
    return
end

currentIndex=cell(length(obj.dividers),1);
for iDim= 1:length(obj.dividers)
    lower=bounds(iDim,1);
    upper=bounds(iDim,2);
    % find those above the lower bound
    E=obj.dividers{iDim}>=lower;
    % find those above the upper bound
    F=obj.dividers{iDim}<=upper;
    % bitshift and add it to the list of those already there
    currentIndex{iDim}=unique([currentIndex{iDim},find(E(2:end)&F(1:(end-1)))]);
   % if isempty(currentIndex{iDim})||E(1)||F(end)
    if isempty(currentIndex{iDim})
        error=1;
    end
end

MX = zeros(1,0);
for iDim=1:length(currentIndex)
    MX=[repmat(MX,length(currentIndex{iDim}),1),reshape(repmat(currentIndex{iDim},size(MX,1),1),[],1)];
end
Multiples=ones(length(currentIndex),1);
for i = 1:(length(currentIndex)-1)
    Multiples(i,1)=prod(obj.nrOfSegments((i+1):end));
end
cells=(MX-1)*Multiples;


if isempty(currentIndex)||error
    cells = 0;
elseif giveAnswerAsIndices
    cells = s2i(obj,MX);
    
    if sum(bounds(:,1)<obj.intervals(:,1))||sum(bounds(:,2)>obj.intervals(:,2))
        cells = [cells,0];
    end
else
    cells=MX;
    if sum(bounds(:,1)<obj.intervals(:,1))||sum(bounds(:,2)>obj.intervals(:,2))
        cells = [cells;zeros(1,length(obj.nrOfSegments))];
    end
end

end

%------------- END OF CODE --------------