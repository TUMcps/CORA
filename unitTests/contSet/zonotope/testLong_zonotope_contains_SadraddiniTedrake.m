function res = testLong_zonotope_contains_SadraddiniTedrake
% testLong_zonotope_contains_SadraddiniTedrake - unit test function
%    of contains using the method of SadraddiniTedrake. We compare the
%    results with an implementation using YALMIP, which is closer to the
%    paper, but implemented less efficiently
%
% Syntax:
%    res = testLong_zonotope_contains_SadraddiniTedrake
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
% See also: -

% Authors:       Matthias Althoff
% Written:       21-July-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% initialize partial results
resPartial = [];

%% create zonotopes -- random cases in 3D (non-zero center)
for iSet = 1:10
    % create zonotope which is supposed to be enclosed (s: small)
    Z_s = zonotope.generateRandom('Dimension',3,'NrGenerators',5);
    
    % create zonotope which is supposed to enclose the smaller one (l: large)
    Z_l = enlarge(zonotope.generateRandom('Dimension',3,'NrGenerators',5), 1.2);
    
    % compute result
    res_original = contains(Z_l, Z_s, 'st');
    
    % compute alternative result
    % Currently, the next function produces an error with the current
    % Yalmip version (as of 09-February-2024), in combination with
    % Matlab>2023b. Thus, instead we compare the results of 'st' with the
    % exact method
    %res_alternative = aux_SadraddiniTedrake_alternative(Z_s, Z_l, 0); % order of Z_l and Z_s is different for this function
    res_alternative = Z_l.contains(Z_s, 'exact');

    if res_alternative && ~res_original
        resPartial(end+1) = false;
    end

    % loop over all types
    %resPartial(end+1) = (res_original == res_alternative);
end

%% create zonotopes -- random cases in 3D (zero center)
for iSet = 1:10
    % create zonotope which is supposed to be enclosed (s: small)
    Z_s = zonotope.generateRandom('Center',zeros(3,1),'Dimension',3,'NrGenerators',5);
    
    % create zonotope which is supposed to enclose the smaller one (l: large)
    Z_l = enlarge(zonotope.generateRandom('Center',zeros(3,1),'Dimension',3,'NrGenerators',5), 1.2);
    
    % compute result
    res_original = contains(Z_l, Z_s, 'st');
    
    % compute alternative result
    % Currently, the next function produces an error with the current
    % Yalmip version (as of 09-February-2024), in combination with
    % Matlab>2023b. Thus, instead we compare the results of 'st' with the
    % exact method
    %res_alternative = aux_SadraddiniTedrake_alternative(Z_s, Z_l, 0); % order of Z_l and Z_s is different for this function
    res_alternative = Z_l.contains(Z_s, 'exact');

    if res_alternative && ~res_original
        resPartial(end+1) = false;
    end
    
    % loop over all types
    %resPartial(end+1) = (res_original == res_alternative);
end

%result of all tests
res = all(resPartial);

end


% Auxiliary functions -----------------------------------------------------

function isIn = aux_SadraddiniTedrake_alternative(Z1, Z2, tol)
% Alternative implementation of the enclosure check by Saddradini and 
% Tedrake

% extract data
x = Z1.center;
y = Z2.center;
X = Z1.generators;
Y = Z2.generators;
nx = size(X, 2);
ny = size(Y, 2);

% yalmip constructs linear programming matrices
Gamma = sdpvar(ny, nx, 'full');
beda = sdpvar(ny, 1);
constraints = [...
    X == Y*Gamma, ...
    y - x == Y*beda, ...
    norm([Gamma, beda], Inf) <= 1+tol];
cost = []; % It suffices to check feasibility here, so no cost function
options = sdpsettings('solver','linprog','verbose',0,'allownonconvex',0); 

% solve linear programming problem
yalmipOptimizer = optimizer(constraints, cost, options, [], {Gamma, beda});

[~, exitFlag] = yalmipOptimizer();

isIn = ~exitFlag;

end

% ------------------------------ END OF CODE ------------------------------
