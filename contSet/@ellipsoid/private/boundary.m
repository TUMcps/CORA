function [Y,L] = boundary(E,N)
% boundary - Deterministically computes N uniformly distributed points on 
%    the boundary of an ellipsoid
%
% Syntax:
%    [Y,L] = boundary(E,N)
%
% Inputs:
%    E - ellipsoid object
%    N - number of samples
%
% Outputs:
%    Y - boundary points
%    L - sampled directions (unit vectors)
%
% Example: 
%    E = ellipsoid([1,0;0,1/2],[1;1]);
%    Y = randPoint(E,1000,'extreme');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ellipsoid/randPoint_

% Authors:       Victor Gassmann
% Written:       13-March-2019
% Last update:   12-March-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if N==0
    throw(CORAerror('CORA:wrongValue','second','positive integer value'));
end

% if only center remains
if rank(E)==0
    Y = repmat(E.q,1,N);
    L = zeros(dim(E),N);
    return;
end

x_rem = [];
[T,~,~] = svd(E.Q);
% shift to zero, then transform to diagonal 
E_nd = T'*(E + (-E.q));

% check if degenerate, and if so, project to lower-dim space
if ~isFullDim(E)
    nt = rank(E_nd);
    x_rem = E_nd.q(nt+1:end);
    E_nd = project(E_nd,1:nt);
end 
n_nd = dim(E_nd);
% compute equally space points on boundary of n_nd-sphere
if n_nd>=2
    L_nd = eq_point_set(n_nd-1,N);
else
    if mod(N,2) == 0
        L_nd = linspace(-1,1,N);
    else
        % ensure that L_nd does not become 0
        L_nd = linspace(-1,1,N+1);
        L_nd = L_nd(1:end-1);
    end
end

Y_nd = zeros(n_nd,N);
% compute boundary point based on support function
for i=1:N
    Y_nd(:,i) = E_nd.Q*L_nd(:,i)/sqrt(L_nd(:,i)'*E_nd.Q*L_nd(:,i));
end
Y_t = [Y_nd;repmat(x_rem,1,N)];
L_t = [L_nd;zeros(length(x_rem),N)];

% backtransform and shift
Y = T*Y_t + E.q;
L = T*L_t;

% ------------------------------ END OF CODE ------------------------------
