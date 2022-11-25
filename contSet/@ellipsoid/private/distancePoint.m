function D = distancePoint(E,Y)
% distancePoint - computes the distance from an ellipsoid to an (array of)
% points
%
% Syntax:  
%    D = distancePoint(E,Y)
%
% Inputs:
%    E - ellipsoid object
%    Y - point or cell array of points
%
% Outputs:
%    D - distance(s) between E and Y
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Victor Gassmann
% Written:      08-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
if iscell(Y)
    Y = cell2mat(reshape(Y,[1,numel(Y)]));
end
n = dim(E);
N = size(Y,2);
x_rem = zeros(0,1);
if E.isdegenerate
    nt = rank(E);
    if nt==0
        D = sqrt(sum((E.q-Y).^2,1));
        return;
    end
    [T,~,~] = svd(E.Q);    
    E = T'*E;
    Y = T'*Y;
    x_rem = E.q(nt+1:end);
    E = project(E,1:nt);
end

% handle nondegenerate case
n_nd = dim(E);
x_nd = sdpvar(n_nd,1);
x = [x_nd;x_rem];
C = (x_nd-E.q)'*inv(E.Q)*(x_nd-E.q) <= 1;
sdpOpts = sdpsettings('verbose',0);
if N>1
    % construct optimizer object to increase computation speed
    y = sdpvar(n,1);
    f_obj = norm(x-y);
    dist = optimizer(C,f_obj,sdpOpts,y,f_obj);
    D = zeros(1,N);
    for i=1:N
        D(i) = dist(Y(:,i));
    end
else
    % only one point, no need for optimizer
    f_obj = norm(x-Y);
    optimize(C,f_obj,sdpOpts);
    % extract optimal point
    D = value(f_obj);
end
%------------- END CODE --------------