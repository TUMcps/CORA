function res = testMP_exponential_Krylov_projected_linSysInput(~)
% testMP_exponential_Krylov_projected_linSysInput - unit_test_function for checking
%    the Krylov method for propagating a zonotope during reachability analysis.
%
% Syntax:
%    res = testMP_exponential_Krylov_projected_linSysInput(~)
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

% enable access to private Krylov functions
currDir = pwd;
cd(strcat(CORAROOT,filesep,"contDynamics",filesep,"@linearSys",filesep,"private"));

% set options --------------------------------------------------------------
options.timeStep = 0.01; %time step size for reachable set computation
options.t = options.timeStep;
options.tFinal = 5;
options.krylovError = 3*eps;
options.krylovStep = 10;

%--------------------------------------------------------------------------

% if 'priv_subspace_Krylov_jaweckiBound works this function is pretty safe
% to be correct, so we only test 5 instances
for i = 1:5
    % create system
    n = randi([30,80],1,1);

    sys = linearSys.generateRandom("StateDimension",n, ...
        'RealInterval',interval(-10,-5),'ImaginaryInterval',interval(-0.5,0.5));
    A = sys.A;
    C = eye(size(A,1));

    % limit number of generators to keep test time short
    Z = zonotope.generateRandom("Dimension",n,"NrGenerators",n);

    c = center(Z);
    G = generators(Z);

    % compute subspace + outer approximation
    [V_c,H_c,V_g,H_g,~] = priv_subspace_Krylov_jaweckiBound(A,Z,options);
    
    % create linSys objects
    c_sys = linearSys(H_c,V_c');
    g_sys = cell(size(G,2));
    for j = 1:size(G,2)
        g_sys{j} = linearSys(H_g{j},V_g{j}');
    end
    
    % call function to propagate initial set via exponential matrix 
    % in Krylov subspace
    [outer_approx,~] = priv_exponential_Krylov_projected_linSysInput(c_sys,g_sys,Z,sys.nrOfDims,options);

    % compute exact solution
    exact_sol = expm(A*options.timeStep)*Z;

    % generate random points at the edges of Z_sol and check if they are in
    % the approximate solution
    p_rand = randPoint(exact_sol,50,'extreme');

    % check point containment
    % for j = size(p_rand,2)
        assertLoop(contains(outer_approx,p_rand,'approx',1e-8),i,j);
    % end
end

% revoke access to private functions
cd(currDir);

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
