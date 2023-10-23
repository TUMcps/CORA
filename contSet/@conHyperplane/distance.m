function val = distance(hyp,S)
% distance - computes the distance from a constrained hyperplane to a set
%    or a point
%
% Syntax:
%    val = distance(hyp,S)
%
% Inputs:
%    obj - conHyperplane object
%    S - contSet object or single point
%
% Outputs:
%    val - distance between constrained hyperplane and set
%
% Example: 
%    hyp = conHyperplane(halfspace([1;1],0),[1 0;-1 0],[2;2])
%    p = [1;1];
%    distance(hyp,p)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Victor Gassmann
% Written:       08-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% compute distance
if isa(S,'double')
    val = distancePoint(hyp,S);
else
    throw(CORAerror('CORA:noops',hyp,S));
end

% ------------------------------ END OF CODE ------------------------------
