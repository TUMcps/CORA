function [res,cert,scaling] = priv_zonotopeContainment_SadraddiniTedrake(Z1, Z2, tol, scalingToggle)
% priv_zonotopeContainment_SadraddiniTedrake - Solves the containment problem using the
%    method described in [1, Corollary 4]. Here, we transform the linear
%    progam into the form required for linprog. A detailed derivation of
%    the transformation can be found in the complementary documentation of
%    CORA.
%
% Syntax:
%    [res,cert,scaling] = priv_zonotopeContainment_SadraddiniTedrake(Z1, Z2, tol, scalingToggle)
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
%           For priv_zonotopeContainment_SadraddiniTedrake, this is an upper bound.
%           Note that computing this scaling factor may significantly
%           increase the runtime.
%
% Example:
%    Z1 = zonotope([0.5 2 3 0;0.5 2 0 3]);
%    Z2 = zonotope([0 -1 1 0; 0 1 0 1]);
%    Z3 = Z2 + [3;0];
%
%    % The function priv_zonotopeContainment_SadraddiniTedrake is called implicitly by
%    % contains
%    contains(Z1,Z2,'approx:st')
%    contains(Z1,Z3,'approx:st')
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
%    [1] Sadraddini et. al: Linear Encodings for Polytope Containment
%        Problems, CDC 2019
%    [2] Kulmburg A., Schäfer L., Althoff M.: 
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

% Extract data; 'circum' stands for 'circumbody', i.e., the zonotope
% that is supposed to contain the 'inbody'
% Extract the generators
G_inbody = Z1.generators;
G_circum = Z2.generators;
% Extract the centers
c_inbody = Z1.center;
c_circum = Z2.center;
% Add the difference of the centers to the generators of the inbody;
% one can show that this yields an equivalent optmization problem
G_inbody = [G_inbody c_inbody-c_circum];

% Extract the number of generators
m_inbody = size(G_inbody,2);
m_circum = size(G_circum,2);

% Extract the dimension
n = size(G_inbody,1);

% Sparsifying
G_inbody = sparse(G_inbody);
G_circum = sparse(G_circum);

% Setting up constraints for the linear program
% Gamma - GammaAux <= 0
A_positive = [speye(m_inbody*m_circum) -speye(m_inbody*m_circum) sparse(m_inbody*m_circum,1)];
b_positive = sparse(m_inbody*m_circum,1);
% -Gamma - GammaAux <= 0
A_negative = [-speye(m_inbody*m_circum) -speye(m_inbody*m_circum) sparse(m_inbody*m_circum,1)];
b_negative = sparse(m_inbody*m_circum,1);
% Row-wise sum of GammaAux <= w
A_rowwiseSum = [sparse(m_circum,m_inbody*m_circum) kron(sparse(ones([1 m_inbody])),speye(m_circum)) -sparse(ones([m_circum 1]))];
b_rowwiseSum = sparse(m_circum,1);

% Building the inequality conditions
A = [A_positive;A_negative;A_rowwiseSum];
b = [b_positive;b_negative;b_rowwiseSum];

% G_circum * Gamma = G_inbody
Aeq = [kron(speye(m_inbody), G_circum) sparse(n*m_inbody,m_inbody*m_circum) sparse(n*m_inbody,1)];
beq = G_inbody(:);

% Defining cost
cost = [sparse(m_inbody*m_circum,1); sparse(m_inbody*m_circum,1); 1];


problem.f = cost;
problem.Aeq = Aeq;
problem.beq = beq;
problem.Aineq = A;
problem.bineq = b;
problem.lb = [];
problem.ub = [];

% solve linear program
[~,scaling,~] = CORAlinprog(problem);

exitflag = scaling <= 1+tol;

if scalingToggle
    res = scaling<=1+tol;
    if res
        cert = true;
    elseif scaling > sqrt(m_inbody)
        % (see Theorem III.4 in [2])
        cert = true;
    else
        cert = false;
    end
else
    res = (exitflag == 1);
    cert = res;
    scaling = NaN;
end
end

% ------------------------------ END OF CODE ------------------------------
