function E = andEllipsoidIA(E,E_c)
% andEllipsoidIA - Computes the inner approximation of the intersection between two
%           ellipsoids
%
% Syntax:  
%    E = andEllipsoidIA(E1,E2)
%
% Inputs:
%    E1 - ellipsoid object
%    E_c - ellipsoid object or array of ellipsoid objects
%
% Outputs:
%    E - ellipsoid after intersection
%
% References:
%   [1] Convex Optimization, Boyd et al.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: andEllipsoidOA

% Author:       Victor Gassmann
% Written:      09-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
E_c = [{E};E_c];
ind_dMask = cellfun(@(e)e.isdegenerate,E_c);
if sum(ind_dMask)>1
    error('At most one ellipsoid can be degenerate!');
end

% make at most last array element degenerate
Ed_c = E_c(ind_dMask);
End_c = E_c(~ind_dMask);
E_c = [End_c;Ed_c];

n = length(E_c{1}.q);
N = length(E_c);
I = eye(n);
T = eye(n);
x_rem = zeros(0,1);
if E_c{end}.isdegenerate
    [T,~,~] = svd(E_c{end}.Q);
    nt = E_c{end}.rank;
    % transform all ellipsoids
    E_c = cellfun(@(e)T'*e,E_c,'Uni',false);
    % project
    x_rem = E_c{end}.q(nt+1:end);
    E_c{end} = project(E_c{end},1:nt);
    % resolve x_rem in E_c{1:end-1} by cutting with hyperplanes 
    % I(nt+1:end,:)*xt = x_rem
    for i=1:N-1
        % intersect each cell element with all hyperplanes
        for j=1:(n-nt)
            Hi = conHyperplane(I(nt+j,:),x_rem(j));
            E_c{i} = E_c{i}&Hi;
            % check if intersection is empty; if yes, overall result is
            % empty
            if isempty(E_c{i})
                E = ellipsoid;
                return;
            end
        end
        % E_c{i} now also has zeros at nt+1:n
        E_c{i} = project(E_c{i},1:nt);
    end
end

n_nd = E_c{1}.dim;
% formalize and solve optimization problem
B = sdpvar(n_nd);
l = sdpvar(N,1);
d = sdpvar(n_nd,1);
f_obj = -geomean(B);
C = [];
for j=1:N
    Aiinv = E_c{j}.Q;
    bi = -inv(E_c{j}.Q)*E_c{j}.q;
    ci = E_c{j}.q'*inv(E_c{j}.Q)*E_c{j}.q - 1;
    % formulate constraint
    Ci = [-l(j)-ci+bi'*Aiinv*bi,zeros(1,n_nd),       (d+Aiinv*bi)';
          zeros(n_nd,1),           l(j)*eye(n_nd),    B;
          d+Aiinv*bi,           B,              Aiinv];
    C = [C,Ci>=0];
end
sol = optimize([C,l>=0],f_obj,sdpsettings('verbose',0));
if sol.problem==0
    Qt = value(B)^2;
    qt = value(d);
    % backtransform
    E = T*ellipsoid([Qt,zeros(n_nd,n-n_nd);zeros(n-n_nd,n)],[qt;x_rem]);
elseif sol.problem==1
    E = ellipsoid;
else
    throw(errOptNumIssue());
end
%------------- END OF CODE --------------