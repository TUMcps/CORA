function I = interval(zB)
% interval - converts a zonotope bundle to an interval according to
%    Proposition 6 in [1]
%
% Syntax:
%    I = interval(zB)
%
% Inputs:
%    zB - zonoBundle object
%
% Outputs:
%    I - interval object
%
% References:
%    [1] M. Althoff. "Zonotope bundles for the efficient computation of 
%        reachable sets", 2011
%
% Other m-files required: interval(constructor)
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       10-November-2010
% Last update:   25-July-2016 (intervalhull replaced by interval)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if zB.parallelSets == 0
    I = interval.empty(dim(zB)); return
end

% enclose all zonotopes by an interval
IHtmp = cell(zB.parallelSets,1);
for i=1:zB.parallelSets
    IHtmp{i} = interval(zB.Z{i});
end

% intersect interval hulls
I = IHtmp{1};
for i=2:zB.parallelSets
    I = and_(I,IHtmp{i},'exact');
end

% ------------------------------ END OF CODE ------------------------------
