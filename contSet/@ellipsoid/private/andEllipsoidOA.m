function E = andEllipsoidOA(E1,E2)
% andEllipsoidOA - Computes the outer approximation of the intersection between two
%           ellipsoids
%
% Syntax:  
%    E = andEllipsoidOA(E1,E2)
%
% Inputs:
%    E1 - ellipsoid object
%    E2 - ellipsoid object or numerical vector
%
% Outputs:
%    E - ellipsoid after intersection
%
% References:
%   [1] Kurzhanskiy, A.A. and Varaiya, P., 2006, December. Ellipsoidal toolbox (ET).
% In Proceedings of the 45th IEEE Conference on Decision and Control (pp. 1498-1503). IEEE.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      09-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
if E1.isdegenerate && E2.isdegenerate
    error('At least one ellipsoid has to be full-dimensional!')
end

%check if ellipsoids are equal (within tolerance)
if E1==E2
    E = E1;
    return;
end
% check if ellipsoids are intersecting
if ~isIntersecting(E1,E2)
    E = ellipsoid;
    return;
end

% make sure that at most E2 is degenerate 
if E1.isdegenerate
    tmp = E2;
    E2 = E1;
    E1 = tmp;
end

n = length(E1.q);
I = eye(n);
T = eye(n);
x2_rem = [];

% handle degeneracy
if E2.isdegenerate
    nt = rank(E2);
    % already handled in ellipsoid/and, but include for safety
    if nt==0
        
    end
    [T,~,~] = svd(E2.Q);
    E2 = T'*E2;
    E1 = T'*E1;
    % project
    x2_rem = E2.q(nt+1:end);
    E2 = project(E2,1:nt);
    % resolve x2_rem in E1 by cutting with hyperplanes 
    % I(nt+1:end,:)*xt = x_rem
    for i=1:(n-nt)
        Hi = conHyperplane(I(nt+i,:)',x2_rem(i));
        E1 = E1 & Hi;
        % since they are intersecting, E1 will not be empty
        assert(~isempty(E1),'Bug: Intersection should not be empty');
    end
    % E1 also has zeros at nt+1:n
    E1 = project(E1,1:nt);
end

n_nd = length(E1.q);
W1 = inv(E1.Q);
W2 = inv(E2.Q);
q1 = E1.q;
q2 = E2.q;
p = compIntersectionParam(W1,q1,W2,q2);

[~,Qt,qt] = rootfnc(p,W1,q1,W2,q2);
if any(eig(Qt)<0)
    E = ellipsoid;
else
    E = T*ellipsoid([Qt,zeros(n_nd,n-n_nd);zeros(n-n_nd,n)],[qt;x2_rem]);
end
%------------- END OF CODE --------------