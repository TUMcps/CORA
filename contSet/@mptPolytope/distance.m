function Val = distance(P,S)
% enclose - Computes the closest distance from set S to Y
%
% Syntax:  
%    P = distance(P,S)
%
% Inputs:
%    P - polytope object
%    S - contSet object
%
% Outputs:
%    Val - distance to S
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/enclose

% Author:       Victor Gassmann
% Written:      26-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
% check dimensions
if ~isa(S,'double')
    if dim(S)~=dim(P)
        error('Dimensions of arguments do not match!');
    end
else
    if size(S,1)~=dim(P)
        error('Dimension of array of points and mptPolytope must match!');
    end
end
    
%% Different distances
if isa(S,'ellipsoid')
    Val = distance(S,P);
elseif isa(S,'double')
    Val = distanceDouble(P,S);
else 
    error('Not implemented yet!');
end
end

%-- helper
function Val = distanceDouble(P,Y)
[n,N] = size(Y);
A = P.P.A;
b = P.P.b;
x = sdpvar(n,1);
C = A*x<=b;
sdpOpts = sdpsettings('verbose',0);
if N>1
    % construct optimizer object to increase computation speed
    y = sdpvar(n,1);
    f_obj = norm(x-y);
    dist = optimizer(C,f_obj,sdpOpts,y,f_obj);
    Val = zeros(1,N);
    for i=1:N
        Val(i) = dist(Y(:,i));
    end
else
    % only one point, no need for optimizer
    f_obj = norm(x-Y);
    optimize(C,f_obj,sdpOpts);
    % extract optimal point
    Val = value(f_obj);
end
end
%------------- END OF CODE --------------