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

% Authors:       Tobias Ladner, Mark Wetzlinger
% Written:       18-September-2023
% Last update:   14-July-2024 (MW, support vertex representation, special case)
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

% copy polytope
P_out = polytope(P);

if N == dim(P)
    % only shuffle dimensions
    P_out = aux_shuffleDimensions(P_out,dims);
    return
end

% halfspace representation (note: we need to check the H representation
% first, since there is a lift-operation which kills the vertices)
if P.isHRep.val

    % lift polytope to higher dimension
    P_out = lift_(P_out,N,dims);
    % override properties
    P_out.bounded.val = P.bounded.val;
    
    % bound new dimensions at 0 using equality constraints
    AeNewDims = eye(N);
    AeNewDims = AeNewDims(~ismember(1:N,dims),:);
    beNewDims = zeros(size(AeNewDims,1),1);
    
    % append to equality constraints
    P_out.Ae_.val = [P_out.Ae_.val; AeNewDims];
    P_out.be_.val = [P_out.be_.val; beNewDims];
end

% vertex representation
if P.isVRep.val
    % add zeros for other dimensions
    V = zeros(N, size(P.V_.val,2));
    V(dims,:) = P.V_.val;
    P_out.V_.val = V;
    P_out.isVRep.val = true;
end

% since we add at least one 'flat' dimension (otherwise we would have
% entered a different path above), the new set is degenerate
P_out.fullDim.val = false;

end


% Auxiliary functions -----------------------------------------------------

function P_out = aux_shuffleDimensions(P_out,dims)
% the 'higher'-dimensional space has the same number of dimensions as P
% only check for the order of the dimensions

if P_out.isHRep.val
    % re-order columns of constraint matrices
    P_out.A_.val = P_out.A_.val(:,dims);
    P_out.Ae_.val = P_out.Ae_.val(:,dims);
end

if P_out.isVRep.val
    % re-order rows of vertices
    P_out.V_.val = P_out.V_.val(dims,:);
end

end

% ------------------------------ END OF CODE ------------------------------
