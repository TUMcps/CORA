function p = gridPoints(I,segments)
% gridPoints - computes uniformly partitioned grid points of an interval;
%
% Syntax:
%    coordinateMat = gridPoints(I,segments)
%
% Inputs:
%    I - interval object
%    segments - number of segments for each dimension (scalar or vector),
%               has to fit dimension of interval!
%
% Outputs:
%    p - (input: nx1 or 1xn interval) matrix with grid points
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

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       22-July-2016
% Last update:   08-November-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input arguments
inputArgsCheck({{I,'att','interval'};
                {segments,'att','numeric',{'nonnan','positive'}}});

% empty case
if representsa_(I,'emptySet',eps)
    p = zeros(dim(I),0);
    return;
end

% check whether interval or interval matrix
[n,m] = size(I);
intTrans = false;
intMat = false;
if n > 1 && m > 1
    intMat = true;
    % make interval matrix a vector for the computation
    I = interval(I.inf(:),I.sup(:));
    nOrig = n; mOrig = m;
    n = n*m;
    segments = segments(:);
elseif m > 1
    % transpose interval if necessary
    intTrans = true;
    I = I';
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
    throw(CORAerror('CORA:wrongValue','second',"larger than 0"));
end
if ~isscalar(segments) && isvector(segments) && length(segments) ~= n
    throw(CORAerror('CORA:wrongValue','second',...
        "a scalar or a vector and fit the dimension of interval"));
end

% radius, center, and number of dimensions with non-zero length
r = rad(I);
c = center(I);

% case: radius of zero-length in all dimensions
if ~any(r)
    if intMat
        % interval matrices... stored in cell-array
        p = {c};
    else
        % regular intervals... stored in array
        if intTrans
            p = c';
        else
            p = c;
        end
    end
    return;
end

% case: non-empty radius in some dimensions
% dims where segments == 1 or radius is zero (same for computation)
segEq1 = segments == 1 | ~any(r,2);

% obtain segment length
segLengthVec = (I.sup-I.inf) ./ (segments-1);
segLengthVec(segEq1) = 0;

% starting point: either at infimum or in the center (if corresponding
% segments-entry is 1)
startingPoint = I.inf;
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
p = zeros(n,nrPerms);
% obtain first grid point
p(:,1) = startingPoint;

% add further grid points
for i = 2:nrPerms
    % shift curr_comb by 1 for correct addition of segLengthVec
    p(:,i) = startingPoint + (perm(i,:)'-1) .* segLengthVec;
end

% transpose if original interval was transposed
if intTrans
    p = p';
end

% switch representation of points back to interval matrix if necessary
if intMat
    temp = cell(nrPerms,1);
    for i = 1:nrPerms
        temp{i} = reshape(p(:,i),nOrig,mOrig);
    end
    p = temp;
end

% ------------------------------ END OF CODE ------------------------------
