function res = testLong_zonotope_minkDiff_RaghuramanKoeln
% testLong_zonotope_minkDiff_RaghuramanKoeln - unit test function of 
%    minkDiff using the method of RaghuramanKoeln. We compare the results
%    with an implementation using YALMIP, which is closer to the paper, but
%    implemented less efficiently
%
% Syntax:
%    res = testLong_zonotope_minkDiff_RaghuramanKoeln
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       27-July-2022
% Last update:   ---
% Last revision: 06-March-2024 (TL, work around yalmip x0 bug)

% ------------------------------ BEGIN CODE -------------------------------

% define small box
smallBox = zonotope([[0;0;0],1e-6*eye(3)]);

% initialize partial results
resvec = true(0);

%% create zonotopes -- random cases in 3D
for iSet = 1:20
    % create minuend
    Z_m = zonotope.generateRandom('Dimension',3,'NrGenerators',5);
    
    % create subtrahend
    Z_s = enlarge(zonotope.generateRandom('Dimension',3,'NrGenerators',5), 0.2);

    if iSet > 10
        % test with zero center
        Z_m = Z_m - Z_m.c;
        Z_s = Z_s - Z_s.c;
    end
    
    % compute result
    Z_res_original = minkDiff(Z_m, Z_s, 'inner:RaghuramanKoeln');
    
    % compute alternative result
    Z_res_alternative = aux_RaghuramanKoeln_alternative(Z_m, Z_s); 
    
    % check whether Minkowski difference returns the empty set
    if representsa(Z_res_original,'emptySet')
        % check if polytope solution is empty as well
        resvec(end+1,1) = representsa(Z_res_alternative,'emptySet');
    elseif representsa(Z_res_alternative,'emptySet')
        % due to a yalmip bug (x0), the alternative solution is not able to
        % compute a solution, and thus aux_RaghuramanKoeln_alternative
        % returns the empty set.
        % We accept that case here...
        resvec(end+1,1) = true;
    else
        % enclosure check: alternative in original + smallBox
        resvec(end+1,1) = contains(Z_res_original + smallBox, Z_res_alternative);

        % enclosure check: alternative in original + smallBox
        resvec(end+1,1) = contains(Z_res_alternative + smallBox, Z_res_original);
    end
    
    if resvec(end) ~= 1
        % MPT toolbox often wrongfully returns that a polytope is empty
        throw(CORAerror('CORA:testFailed'))
    end
end

%result of all tests
res = all(resvec);

end


% Auxiliary functions -----------------------------------------------------

function Z = aux_RaghuramanKoeln_alternative(Z1, Z2)
% Computes the Minkowski difference according to [2]
% Z1: minuend
% Z2: subtrahend
   
    % extract data
    c1 = Z1.c;
    c2 = Z2.c;
    G1 = Z1.G;
    G2 = Z2.G;
    gens1 = size(G1, 2); % number of generators of Z1
    gens2 = size(G2, 2); % number of generators of Z2
    n = length(c1);

    %% yalmip constructs linear programming matrices
    % instantiate variables
    Gamma = sdpvar(gens1, gens1 + 2*gens2, 'full');
    beta = sdpvar(gens1, 1);
    phi = sdpvar(gens1 + gens2, 1); % diagonal elements of Phi
    cd = sdpvar(n, 1); % shift of center to be optimized
    
    % constraints
    constraints = [...
        [[G1, G2]*diag(phi), G2] == G1*Gamma, ... % eq. (32a) in [2]
        c1 - (cd + c2) == G1*beta, ... % eq. (32b) in [2]
        norm(Gamma, Inf)*ones(gens1 + 2*gens2,1) + norm(beta, Inf) <= 1]; % eq. (32c) in [2]
    
    % objective
    objective = -sum(phi);
    
    options = sdpsettings('solver','linprog', 'verbose',0, 'allownonconvex',0); 

    % solve linear programming problem
    diagnostics = optimize(constraints, objective, options);
    
    % if solution exists
    if diagnostics.problem == 0
        % optimal solution of phi
        phi_opt = value(phi);
        cd_opt = value(cd);

        % result
        Z = zonotope([cd_opt,[G1, G2]*diag(phi_opt)]);
    else
    % no solution exists
        Z = zonotope.empty(n);
    end
end

% ------------------------------ END OF CODE ------------------------------
