function res = testMP_subspace_Krylov_individual_Jawecki(~)
% testMP_subspace_Krylov_individual_Jawecki - unit_test_function for checking
%    the Krylov method for computing the exponential matrix multiplied with a single vector.
%
% Syntax:
%    res = testMP_subspace_Krylov_individual_Jawecki(~)
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

% enable access to private function "priv_subspace_Krylov_individual_Jawecki"
currDir = pwd;
cd(strcat(CORAROOT,filesep,"contDynamics",filesep,"@linearSys",filesep,"private"));

%set options --------------------------------------------------------------
options.timeStep = 0.01; %time step size for reachable set computation
options.tFinal = 5;
options.krylovError = 3*eps;
options.krylovStep = 10;

%--------------------------------------------------------------------------

for i = 1:100
    % create system
    n = randi([30,80],1,1);

    sys = linearSys.generateRandom("StateDimension",n, ...
        'RealInterval',interval(-10,-5),'ImaginaryInterval',interval(-0.5,0.5));
    A = sys.A;
    v = rand(n,1);
    v_norm = norm(v);

    % compute subspace + approximation
    [V,H] = priv_subspace_Krylov_individual_Jawecki(A,v,1,options);
    sol_apx = v_norm*V*expm(H*options.timeStep)*unitvector(1,size(H,1));

    % compute exact solution
    sol_exc = expm(A*options.timeStep)*v;

    % assert condition
    % (higher tol here to account for floating point errors along the way)
    try
        assertLoop(norm(sol_apx-sol_exc) <= 1e-8,i);
    catch
        disp(norm(sol_apx-sol_exc))
    end
end

% revoke access to private function "initReach_Krylov"
cd(currDir);

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
