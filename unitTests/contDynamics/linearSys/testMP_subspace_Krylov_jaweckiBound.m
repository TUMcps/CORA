function res = testMP_subspace_Krylov_jaweckiBound(~)
% testMP_subspace_Krylov_jaweckiBound - unit_test_function for checking
%    the Krylov method for computing the exponential matrix multiplied with a zonotope.
%
% Syntax:
%    res = testMP_subspace_Krylov_jaweckiBound(~)
%
% Inputs:
%    no
%
% Outputs:
%    res - true/false

% Authors:       Maximilian Perschl
% Written:       28-May-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
% enable access to private functions "priv_subspace_Krylov_individual_Jawecki"
% enable access to private functions "priv_subspace_Krylov_jaweckiBound"
currDir = pwd;
cd(strcat(CORAROOT,filesep,"contDynamics",filesep,"@linearSys",filesep,"private"));

%set options --------------------------------------------------------------
options.timeStep = 0.01; %time step size for reachable set computation
options.tFinal = 5;
options.krylovError = 3*eps;
options.krylovStep = 10;

%--------------------------------------------------------------------------

for i = 1:20
    % create system
    n = randi([30,80],1,1);

    sys = linearSys.generateRandom("StateDimension",n, ...
        'RealInterval',interval(-10,-5),'ImaginaryInterval',interval(-0.5,0.5));
    A = sys.A;
    
    % limit number of generators to keep test time short
    Z = zonotope.generateRandom("Dimension",n,"NrGenerators",n);

    c = center(Z);
    G = generators(Z);

    % compute subspace + outer approximation
    [V_c,H_c,V_g,H_g,~] = priv_subspace_Krylov_jaweckiBound(A,Z,options);
    
    % compute outer approximate solution set
    c_new = norm(c)*V_c*expm(H_c*options.timeStep)*unitvector(1,size(H_c,1));
    G_new = zeros(size(G));
    for j = 1:size(G,2)
        g = G(:,j);
        G_new(:,j) = norm(g)*V_g{j}*expm(H_g{j}*options.timeStep)*unitvector(1,size(H_g{j},1));
    end
    Z_approx = zonotope(c_new,G_new);

    % compute exact solution
    Z_sol = expm(A*options.timeStep)*Z;

    % generate random points at the edges of Z_sol and check if they are in
    % the approximate solution
    p_rand = randPoint(Z_sol,50,'extreme');

    % check point containment
    % for j = size(p_rand,2)
        assertLoop(contains(Z_approx,p_rand,'approx',1e-8),i,j);
    % end
end

% revoke access to private function "initReach_Krylov"
cd(currDir);

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
