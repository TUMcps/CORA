function coordinateMat = gridPoints(obj,segments)
% gridPoints - computes uniformly partitioned grid points of an interval;
%
% Syntax:  
%    coordinateMat = gridPoints(obj,segments)
%
% Inputs:
%    obj - interval object
%    segments - number of segments for each dimension (scalar or vector)
%
% Outputs:
%    coordinateMat - (input: nx1 or 1xn interval) matrix with grid points
%    	(input: nxm interval) cell array with grid points as matrices
%
% Example: 
%    I = interval([-1;1],[-1;2]);
%    coordinates = gridPoints(I,10);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff, Mark Wetzlinger
% Written:      22-July-2016
% Last update:  08-November-2021
% Last revision:---

%------------- BEGIN CODE --------------

% empty case
if isempty(obj)
    coordinateMat = [];
    return;
end

% check whether interval or interval matrix
[n,m] = size(obj);
intTrans = false;
intMat = false;
if n > 1 && m > 1
    intMat = true;
    % make interval matrix a vector for the computation
    obj = interval(obj.inf(:),obj.sup(:));
    nOrig = n; mOrig = m;
    n = n*m;
    segments = segments(:);
elseif m > 1
    % transpose interval if necessary
    intTrans = true;
    obj = obj';
    n = m;
    segments = segments';
end
% make segments a vector
if isscalar(segments)
    segments = repelem(segments,n)';
end

% check if segments correct
if ~all(segments > 0)
    % segments have to be larger than 0
    [msg,id] = errWrongInput('segments'); error(id,msg);
end
if ~isscalar(segments) && isvector(segments) && length(segments) ~= n
    [msg,id] = errWrongInput('segments'); error(id,msg);
end

% radius, center, and number of dimensions with non-zero length
r = rad(obj);
c = center(obj);

% case: radius of zero-length in all dimensions
if ~any(r)
    if intMat
        % interval matrices... stored in cell-array
        coordinateMat = {c};
    else
        % regular intervals... stored in array
        if intTrans
            coordinateMat = c';
        else
            coordinateMat = c;
        end
    end
    return;
end

% case: non-empty radius in some dimensions
% dims where segments == 1 or radius is zero (same for computation)
segEq1 = segments == 1 | ~any(r,2);

% obtain segment length
segLengthVec = (obj.sup-obj.inf) ./ (segments-1);
segLengthVec(segEq1) = 0;

% starting point: either at infimum or in the center (if corresponding
% segments-entry is 1)
startingPoint = obj.inf;
startingPoint(segEq1) = c(segEq1);
segments(segEq1) = 1;

% obtain segment permutations: only take dimensions with non-zero radius
perm = combinator(max(segments),length(segments),'p','r');
nrPerms = length(perm(:,1));
% note: full_fact_mod does not work for length(segments) > 2
% segments = segments(~segEq1);
% perm = full_fact_mod(segments);

% remove all permutations which exceed segments (only necessary if segments
% are different for different dimensions)
if min(segments) ~= max(segments)
    keepIdx = true(nrPerms,1);
    for i=2:nrPerms
        if any(perm(i,:)' > segments)
            keepIdx(i) = false;
        end
    end
    perm = perm(keepIdx,:);
    nrPerms = length(perm(:,1));
end

% initialize output matrix
coordinateMat = zeros(n,nrPerms);
% obtain first grid point
coordinateMat(:,1) = startingPoint;

% add further grid points
for i = 2:nrPerms
    % shift curr_comb by 1 for correct addition of segLengthVec
    coordinateMat(:,i) = startingPoint + (perm(i,:)'-1) .* segLengthVec;
end

% transpose if original interval was transposed
if intTrans
    coordinateMat = coordinateMat';
end

% switch representation of points back to interval matrix if necessary
if intMat
    temp = cell(nrPerms,1);
    for i = 1:nrPerms
        temp{i} = reshape(coordinateMat(:,i),nOrig,mOrig);
    end
    coordinateMat = temp;
end

%------------- END OF CODE --------------