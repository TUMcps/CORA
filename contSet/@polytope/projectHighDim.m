function P_out = projectHighDim(P,N,dims)
% projectHighDim - lifts a polytope object to a higher-dimensional space, 
%    the added dimensions are bounded with 0
%
% Syntax:
%    P_out = projectHighDim(P,N,dims)
%
% Inputs:
%    P - polytope object
%    N - dimension of the higher-dimensional space
%    dims - states of the higher-dimensional space that correspond to the
%           states of the lower-dimensional polytope object
%
% Outputs:
%    P_out - polytope object in the high-dimensional space
%
% Example: 
%    A = [1 0;-1 0;0 1;0 -1;1 1];
%    b = [1;1;1;1;1];
%    P = polytope(A,b);
%    P_ = projectHighDim(P,10,[4,5]);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polytope/lift

% Authors:       Tobias Ladner
% Written:       18-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
if N < dim(P)
    throw(CORAerror('CORA:wrongValue','second',...
        ['Dimension of higher-dimensional space must be larger than ' ...
        'or equal to the dimension of the given polytope.']));
elseif length(dims) ~= dim(P)
    throw(CORAerror('CORA:wrongValue','third',...
        ['Number of dimensions in higher-dimensional space must match '...
        'the dimension of the given polytope.']));
end

% lift polytope to higher dimension
P_lift = lift(P,N,dims);

% bound new dimensions at 0
AeNewDims = eye(N);
AeNewDims = AeNewDims(~ismember(1:N,dims),:);
beNewDims = zeros(size(AeNewDims,1),1);

% append to equality constraints
P_out = P_lift;
P_out.Ae = [P_out.Ae;AeNewDims];
P_out.be = [P_lift.be;beNewDims];

% inherit unchanged properties
P_out.emptySet.val = P.emptySet.val;
if isempty(P.fullDim.val)
    P_out.fullDim.val = [];
else
    P_out.fullDim.val = P.fullDim.val && N == length(dims);
end
P_out.minHRep.val = P.minHRep.val;
P_out.minVRep.val = P.minVRep.val;
P_out.bounded.val = P.bounded.val;

% add zeros to vertices
if ~isempty(P.V.val)
    V = zeros(N, size(P.V.val,2));
    V(dims,:) = P.V.val;
    P_out.V.val = V;
end

% ------------------------------ END OF CODE ------------------------------
