function res = inEllipsoid(E1,E2)
% inEllipsoid - checks whether ellipsoid E2 is contained in ellipsoid E1
%
% Syntax:  
%    E = inEllipsoid(E1,E2)
%
% Inputs:
%    E1,E2      - ellipsoid object
%
% Outputs:
%    res - true if contained, false otherwise
%
% References:
%            [1] Yildirim, E.A., 2006. On the minimum volume covering 
%                ellipsoid of ellipsoids. SIAM Journal on Optimization, 
%                17(3), pp.621-641.     
%            [2] SDPT3: url: http://www.math.nus.edu.sg/~mattohkc/sdpt3.html
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      16-March-2021
% Last update:  27-July-2021 (allowed both E1,E2 to be degenerate)
% Last revision:---

%------------- BEGIN CODE --------------
if dim(E1) ~= dim(E1)
    error('dimensions do not match');
end
n = dim(E1);

% normalize to prevent numerical issues
% first, get sense of scale for both
max_val = max([svd(E1.Q);svd(E2.Q)]);
E1 = 1/max_val*E1;
E2 = 1/max_val*E2;
q1 = E1.q;
E1 = E1 + (-q1);
E2 = E2 + (-q1);

% if both are degenerate
if E1.isdegenerate && E2.isdegenerate
    nt1 = rank(E1);
    nt2 = rank(E2);
    % cannot be contained if one has more actual dimensions
    if nt1<nt2
        res = false;
        return;
    end
    Q_sum = E1.Q+E2.Q;
    n_nd = rank(ellipsoid(Q_sum));
    % can only be contained if degeneracies are aligned
    if n_nd<n
        [T,~,~] = svd(Q_sum);
        E1t = T'*E1;
        E2t = T'*E2;
        if all(withinTol(E1t.q(n_nd+1:end),E2t.q(n_nd+1:end),E1t.TOL))
            res = inEllipsoid(project(E1t,1:n_nd),project(E2t,1:n_nd));
        else
            res = false;
        end
        return;
    else
        res = false;
        return;
    end
    
end

% if E1 is the degenerate one, E2 cannot be contained
if E1.isdegenerate
    res = false;
    return;
end

TOL = min(E1.TOL,E2.TOL);

if E2.isdegenerate
    nt = rank(E2);
    if nt==0
        % E2 contains only single point
        res = in(E1,E2.q);
        return;
    end
    [T,~,~] = svd(E2.Q);
    E2 = T'*E2;
    E1 = T'*E1;
    % project
    x2_rem = E2.q(nt+1:end);
    E2 = project(E2,1:nt);
    % resolve x2_rem in E1t by cutting with hyperplanes 
    % I(nt+1:end,:)*xt = x_rem
    n = E1.dim;
    I = eye(n);
    for i=1:(n-nt)
        Hi = conHyperplane(I(nt+i,:),x2_rem(i));
        E1 = E1 & Hi;
        % if intersection is empty, E2 cannot be contained in E1
        if isempty(E1)
            res = false;
            return;
        end
    end
    % now, E1 also has "0"s at nt+1:end
    E1 = project(E1,1:nt);
end

% simulatenous diagonalization: Find Tb such that
% Tb'*Q1*Tb = I and Tb'*Q2*Tb = D (diagonal)
% if max(diag(D))<=1 => contained
Q_m12 = inv(sqrtm(E1.Q));
[O,~,~] = svd(Q_m12*E2.Q*Q_m12);
Tb = Q_m12*O';
if all(withinTol(E1.q,E2.q,TOL)) && max(diag(Tb'*E2.Q*Tb))<=1+TOL
    res = true;
    return;
end

Qs = E1.Q;
Q = E2.Q;
qs = E1.q;
q = E2.q;

%For more details, see [1]
t = sdpvar;
Qs_ = inv(Qs);
Q_ = inv(Q);
M = [Q_, -Q_*q;
      -(Q_*q)', q'*Q_*q-1];
Ms = [Qs_, -Qs_*qs;
      -(Qs_*qs)', qs'*Qs_*qs-1];

options = sdpsettings('verbose',0);
%when solving problems with higher-dimensional ellipsoids, sdpt3 [2] is
%recommended to avoid numerical issues
if exist('sdpt3','file')
    options.solver = 'sdpt3';
end
diagnostics = optimize(t*M>=Ms,[],options);
%either feasible or not feasible
if ~any(diagnostics.problem==[1,0])
    throw(errOptNumIssue());
end
if value(t)==0
    error('problem with solution of feasibility problem (t)');
end
res = ~diagnostics.problem;
%------------- END OF CODE --------------