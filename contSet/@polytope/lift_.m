function P_out = lift_(P,N,dims)
% lift_ - lifts a polytope object to a higher-dimensional space, the added
%    dimensions are unbounded
%
% Syntax:
%    P_out = lift_(P,N,dims)
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
%    P_ = lift(P,10,[4,5]);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/lift, polytope/projectHighDim

% Authors:       Niklas Kochdumper
% Written:       06-July-2018
% Last update:   31-July-2023 (MW, support equality constraints)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
if N < dim(P)
    throw(CORAerror('CORA:wrongValue','second',...
        ['Dimension of higher-dimensional space must be larger than '...
        'or equal to the dimension of the given polytope.']));
elseif length(dims) ~= dim(P)
    throw(CORAerror('CORA:wrongValue','third',...
        ['Number of dimensions in higher-dimensional space must match ' ...
        'the dimension of the given polytope.']));
end

% project constraints to higher-dimensional space
A = zeros(size(P.A,1),N);
A(:,dims) = P.A;
Ae = zeros(size(P.Ae,1),N);
Ae(:,dims) = P.Ae;

% construct the resulting high dimensional polytope object
P_out = polytope(A,P.b,Ae,P.be);

% inherit unchanged properties
P_out.emptySet.val = P.emptySet.val;
P_out.fullDim.val = P.fullDim.val;
P_out.minHRep.val = P.minHRep.val;

% check bounded and vertices
if N > length(dims)
    % set cannot be bounded since added dimensions are unbounded
    P_out.bounded.val = false;
    % delete vertices as set is unbounded
    P_out.V.val = [];
    P_out.minVRep.val = [];
else % N == length(dims)
    % not actually lifted
    P_out.bounded.val = P.bounded.val;
    % just shuffle vertices
    if ~isempty(P.V.val)
        P_out.V.val = P.V.val(dims,:);
    end
    P_out.minVRep.val = P.minVRep.val;
end

end

% ------------------------------ END OF CODE ------------------------------
