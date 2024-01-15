function E = andEllipsoidOA(E1,E2)
% andEllipsoidOA - Computes an outer-approximation of the intersection
%    of an ellipsoid and another ellipsoid or a vector
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
%   [1] Kurzhanskiy, A.A. and Varaiya, P., 2006, December. Ellipsoidal
%       toolbox (ET). In Proceedings of the 45th IEEE Conference on
%       Decision and Control (pp. 1498-1503). IEEE.
%   [2] https://www2.eecs.berkeley.edu/Pubs/TechRpts/2006/EECS-2006-46.pdf
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Victor Gassmann
% Written:       09-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if ~isFullDim(E1) && ~isFullDim(E2)
    throw(CORAerror('CORA:degenerateSet',...
        'At least one ellipsoid has to be full-dimensional!'));
end

% check if ellipsoids are equal (within tolerance)
if E1==E2
    E = E1;
    return;
end

% check if ellipsoids are intersecting
if ~isIntersecting_(E1,E2,'exact')
    E = ellipsoid.empty(dim(E1));
    return;
end

% make sure that at most E2 is degenerate 
if ~isFullDim(E1)
    tmp = E2;
    E2 = E1;
    E1 = tmp;
end

% define values if E2 is not degenerate (reused later)
n = length(E1.q);
I = eye(n);
T = eye(n);
x2_rem = [];

% handle degeneracy
if ~isFullDim(E2)
    nt = rank(E2);
    [T,~,~] = svd(E2.Q);
    % transform E1 and E2 into new space where E2 degeneracies are exposed
    % as last n-nt dimensions
    E2 = T'*E2;
    E1 = T'*E1;
    % project
    x2_rem = E2.q(nt+1:end);
    E2 = project(E2,1:nt);
    % resolve x2_rem in E1 by cutting with hyperplanes 
    % I(nt+1:end,:)*xt = x_rem
    for i=1:(n-nt)
        Hi = conHyperplane(I(nt+i,:)',x2_rem(i));
        E1 = and_(E1,Hi,'outer');
        % since they are intersecting, E1 will not be empty
        assert(~representsa_(E1,'emptySet',eps),'Bug: Intersection should not be empty');
    end

    if nt==1
        % 1D case, handle separately
        % in that case, the exact result is E1
        E = T*E1;
        return;
    end

    % E1 also has zeros at nt+1:n
    E1 = project(E1,1:nt);
end

n_nd = length(E1.q);

% define inverted matrices for computation of intersection parameter
W1 = inv(E1.Q);
W2 = inv(E2.Q);
q1 = E1.q;
q2 = E2.q;
% compute intersection parameter (see [2, Sec. 2.2.5])
p = compIntersectionParam(W1,q1,W2,q2);

% find the shape and center of the parameterized family of shapes and
% centers by using the calculated p
[~,Qt,qt] = rootfnc(p,W1,q1,W2,q2);
% check if result describes an ellipsoid; if not => empty set
if any(eig(Qt)<0)
    E = ellipsoid;
else
    E = T*ellipsoid([Qt,zeros(n_nd,n-n_nd);zeros(n-n_nd,n)],[qt;x2_rem]);
end

% ------------------------------ END OF CODE ------------------------------
