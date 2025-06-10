function buckets = getMergeBuckets(obj)
% getMergeBuckets - computes points in the output space of the domain which
%    are suitable for merging similar neurons
%
% Syntax:
%    buckets = getMergeBuckets(obj)
%
% Inputs:
%    obj - nnActivationLayer
%
% Outputs:
%    buckets - numeric
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nnLayer

% Authors:       Tobias Ladner
% Written:       23-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% usually given by extreme points
buckets = obj.f([-inf, inf]);
% filter infinity
buckets = buckets(~isinf(buckets));

% ------------------------------ END OF CODE ------------------------------
