function [res,cert,scaling] = priv_zonotopeContainment_SadraddiniTedrakeDual(Z1, Z2, tol, scalingToggle)
% priv_zonotopeContainment_SadraddiniTedrakeDual - Solves the containment problem using the
%    method described in [1, Theorem 3.5] (dual formulation of
%    priv_zonotopeContainment_SadraddiniTedrake).
%
% Syntax:
%    [res,cert,scaling] = priv_zonotopeContainment_SadraddiniTedrakeDual(Z1, Z2, tol, scalingToggle)
%
% Inputs:
%    Z1 - zonotope object, inbody
%    Z2 - zonotope object, circumbody
%    tol - tolerance for the containment check; the higher the
%       tolerance, the more likely it is that points near the boundary of
%       Z2 will be detected as lying in Z2, which can be useful to
%       counteract errors originating from floating point errors.
%    scalingToggle - if set to 'true', scaling will be computed (see
%       below).
%
% Outputs:
%    res - true/false
%    cert - returns true iff the result of res could be
%           verified. For example, if res=false and cert=true, Z1 is
%           guaranteed to not be contained in Z2, whereas if res=false and
%           cert=false, nothing can be deduced (Z1 could still be
%           contained in Z2).
%           If res=true, then cert=true.
%    scaling - returns the smallest number 'scaling', such that
%           scaling*(Z2 - center(Z2)) + center(Z2) contains Z1.
%           For priv_zonotopeContainment_SadraddiniTedrakeDual, this is an upper bound.
%           Note that computing this scaling factor may significantly
%           increase the runtime.
%
% Example:
%    Z1 = zonotope([0.5 2 3 0;0.5 2 0 3]);
%    Z2 = zonotope([0 -1 1 0; 0 1 0 1]);
%    Z3 = Z2 + [3;0];
%
%    % The function priv_zonotopeContainment_SadraddiniTedrakeDual is called implicitly by
%    % contains
%    contains(Z1,Z2,'approx:stDual')
%    contains(Z1,Z3,'approx:stDual')
%
%    figure; hold on;
%    plot(Z1,[1,2],'b');
%    plot(Z2,[1,2],'g');
%
%    figure; hold on;
%    plot(Z1,[1,2],'b');
%    plot(Z3,[1,2],'r');
%
%
% References:
%    [1] A. Kulmburg, M. Althoff.: Approximability of the Containment
%       Problem for Zonotopes and Ellipsotopes, 2024
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/contains_

% Authors:       Adrian Kulmburg
% Written:       05-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

G_inbody = Z1.generators;
G_circum = Z2.generators;

c_inbody = Z1.center;
c_circum = Z2.center;

G_inbody = [G_inbody c_inbody-c_circum];

m_inbody = size(G_inbody,2);
m_circum = size(G_circum,2);

n = size(G_inbody,1);

% Sparsifying
G_inbody = sparse(G_inbody);
G_circum = sparse(G_circum);

positive_constraint = [kron(speye(m_inbody), G_circum') kron(sparse(ones([m_inbody 1])), -speye(m_circum))];
negative_constraint = [kron(speye(m_inbody), -G_circum') kron(sparse(ones([m_inbody 1])), -speye(m_circum))];
summation = [sparse(1,n * m_inbody) sparse(ones([1 m_circum]))];

A = [positive_constraint;negative_constraint];
b = [sparse(2*m_circum*m_inbody,1)];

Aeq = summation;
beq = 1;

cost = [-G_inbody(:);sparse(m_circum,1)];

problem.f = cost;
problem.Aeq = Aeq;
problem.beq = beq;
problem.Aineq = A;
problem.bineq = b;
problem.lb = [];
problem.ub = [];

% solve linear program
[X,scaling,exitflag] = CORAlinprog(problem);

if exitflag == -3
    % If the problem is unbounded, we need to manually set the scaling
    % to avoid potential errors
    scaling = inf;
    res = false;
    cert = true;
else
    X = reshape(X(1:n*m_inbody), [n m_inbody]);

    scaling = -scaling;
    res = scaling<=1+tol;
    if res
        cert = true;
    elseif scaling > sqrt(m_inbody)
        % (see Theorem III.4 in [1])
        cert = true;
    elseif all(all(ismember(X,X(:,1)) | ismember(X,-X(:,1))))
        % (see Theorem III.5 in [1])
        cert = true;
    else
        cert = false;
    end
end
end

% ------------------------------ END OF CODE ------------------------------
