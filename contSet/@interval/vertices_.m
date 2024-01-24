function V = vertices_(I,varargin)
% vertices_ - Computes vertices of an interval object
%
% Syntax:
%    V = vertices_(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    V - vertices object
%
% Example: 
%    I = interval([1; -1], [2; 1]);
%    V = vertices(I);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/vertices, zonotope/vertices_

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       24-July-2006 
% Last update:   27-March-2023 (MW, rename vertices_)
%                28-April-2023 (MW, remove duplicates in degenerate case)
% Last revision: 05-April-2023 (MW, rewrite to support unbounded intervals)

% ------------------------------ BEGIN CODE -------------------------------

% empty case
if dim(I) == 0
    V = []; return
end

% check whether there is a non-zero radius in all dimensions
idxZeroDim = false(dim(I),1);
if ~isFullDim(I)
    % remove dimensions with zero radius -> save indices and add later
    idxZeroDim = withinTol(rad(I),0);
    valZeroDim = I.inf(idxZeroDim);
    I = project(I,~idxZeroDim);
end

% compute all possible combinations of lower/upper bounds
fac = logical(combinator(2,dim(I),'p','r')-1);
nrComb = size(fac,1);

% init all points with lower bound
V = repmat(I.inf,1,nrComb);
% read out supremum
ub = I.sup;

% loop over all factors
for i=1:nrComb
    V(fac(i,:)',i) = ub(fac(i,:));
end

% add back removed dimensions
if any(idxZeroDim)
    V_ = V;
    V = zeros(dim(I),size(V_,2));
    V(idxZeroDim,:) = valZeroDim;
    V(~idxZeroDim,:) = V_;
end

% ------------------------------ END OF CODE ------------------------------
