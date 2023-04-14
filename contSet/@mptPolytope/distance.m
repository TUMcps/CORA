function val = distance(P,S)
% enclose - Computes the closest distance from set S to Y
%
% Syntax:  
%    P = distance(P,S)
%
% Inputs:
%    P - mptPolytope object
%    S - contSet object
%
% Outputs:
%    val - distance to S
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
            throw(CORAerror('CORA:dimensionMismatch',P,S));
        end
    else
        if size(S,1)~=dim(P)
            throw(CORAerror('CORA:dimensionMismatch',P,S));
        end
    end
        
    % different distances
    if isa(S,'ellipsoid')
        val = distance(S,P);
    elseif isa(S,'double')
        val = distanceDouble(P,S);
    else 
        throw(CORAerror('CORA:noops',P,S));
    end

end


% Auxiliary Function ------------------------------------------------------
function val = distanceDouble(P,Y)

    % read out size of points 
    [n,N] = size(Y);

    % read out halfspaces from polytope
    A = P.P.A;
    b = P.P.b;

    x = sdpvar(n,1);
    C = A*x <= b;
    sdpOpts = sdpsettings('verbose',0);

    % branch depending on number of points
    if N>1
        % construct optimizer object to increase computation speed
        y = sdpvar(n,1);
        f_obj = norm(x-y);
        dist = optimizer(C,f_obj,sdpOpts,y,f_obj);
        val = zeros(1,N);
        for i=1:N
            val(i) = dist(Y(:,i));
        end
    else
        % only one point, no need for optimizer
        f_obj = norm(x-Y);
        optimize(C,f_obj,sdpOpts);
        % extract optimal point
        val = value(f_obj);
    end

end

%------------- END OF CODE --------------