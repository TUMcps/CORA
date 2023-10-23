function res = testLongDuration_zonotope_minkDiff_RaghuramanKoeln
% testLongDuration_zonotope_minkDiff_RaghuramanKoeln - unit test function of 
% minus using the method of RaghuramanKoeln. We compare the results with an
% implementation using YALMIP, which is closer to the paper, but
% implemented less efficiently
%
% Syntax:
%    res = testLongDuration_zonotope_minkDiff_RaghuramanKoeln
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Matthias Althoff
% Written:       27-July-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% define small box
smallBox = zonotope([[0;0;0],1e-6*eye(3)]);

% initialize partial results
resPartial = [];

%% create zonotopes -- random cases in 3D (non-zero center)
for iSet = 1:10
    % create minuend
    Z_m = zonotope.generateRandom('Dimension',3,'NrGenerators',5);
    
    % create subtrahend
    Z_s = enlarge(zonotope.generateRandom('Dimension',3,'NrGenerators',5), 0.2);
    
    % compute result
    Z_res_original = minus(Z_m, Z_s, 'RaghuramanKoeln');
    
    % compute alternative result
    Z_res_alternative = aux_RaghuramanKoeln_alternative(Z_m, Z_s); 
    
    % check whether Minkowski difference returns the empty set
    if isempty(Z_res_original)
        % check if polytope solution is empty as well
        resPartial(end+1) = isempty(Z_res_alternative);
    else
        % enclosure check: alternative in original + smallBox
        resPartial(end+1) = in(Z_res_original + smallBox, Z_res_alternative);

        % enclosure check: alternative in original + smallBox
        resPartial(end+1) = in(Z_res_alternative + smallBox, Z_res_original);
    end
    
    if resPartial(end) ~= 1
        % MPT toolbox often wrongfully returns that a polytope is empty
        path = pathFailedTests(mfilename());
        save(path,'Z_m','Z_s');
        disp(['error in Minkowski difference']);
    end
end

%% create zonotopes -- random cases in 3D (zero center)
for iSet = 1:10
    % create minuend
    Z_m = zonotope.generateRandom('Center',zeros(3,1),'Dimension',3,'NrGenerators',5);
    
    % create subtrahend
    Z_s = enlarge(zonotope.generateRandom('Center',zeros(3,1),'Dimension',3,'NrGenerators',5), 0.2);
    
    % compute result
    Z_res_original = minus(Z_m, Z_s, 'RaghuramanKoeln');
    
    % compute alternative result
    Z_res_alternative = aux_RaghuramanKoeln_alternative(Z_m, Z_s); 
    
    % check whether Minkowski difference returns the empty set
    if isempty(Z_res_original)
        % check if polytope solution is empty as well
        resPartial(end+1) = isempty(Z_res_alternative);
    else
        % enclosure check: alternative in original + smallBox
        resPartial(end+1) = in(Z_res_original + smallBox, Z_res_alternative);

        % enclosure check: alternative in original + smallBox
        resPartial(end+1) = in(Z_res_alternative + smallBox, Z_res_original);
    end
    
    if resPartial(end) ~= 1
        % MPT toolbox often wrongfully returns that a polytope is empty
        path = pathFailedTests(mfilename());
        save(path,'Z_m','Z_s');
        disp(['error in Minkowski difference']);
    end
end

%result of all tests
res = all(resPartial);

end


% Auxiliary functions -----------------------------------------------------

function Z = aux_RaghuramanKoeln_alternative(Z1, Z2)
% Computes the Minkowski difference according to [2]
% Z1: minuend
% Z2: subtrahend
   
    % extract data
    c1 = Z1.center;
    c2 = Z2.center;
    G1 = Z1.generators;
    G2 = Z2.generators;
    gens1 = size(G1, 2); % number of generators of Z1
    gens2 = size(G2, 2); % number of generators of Z2
    dim = length(c1);

    %% yalmip constructs linear programming matrices
    % instantiate variables
    Gamma = sdpvar(gens1, gens1 + 2*gens2, 'full');
    beta = sdpvar(gens1, 1);
    phi = sdpvar(gens1 + gens2, 1); % diagonal elements of Phi
    cd = sdpvar(dim, 1); % shift of center to be optimized
    
    % constraints
    constraints = [...
        [[G1, G2]*diag(phi), G2] == G1*Gamma, ... % eq. (32a) in [2]
        c1 - (cd + c2) == G1*beta, ... % eq. (32b) in [2]
        norm([Gamma, beta], Inf) <= 1]; % eq. (32c) in [2]
    
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
        Z = [];
    end
end

% ------------------------------ END OF CODE ------------------------------
