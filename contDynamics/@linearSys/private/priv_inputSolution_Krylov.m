function linsys = priv_inputSolution_Krylov(linsys,params,options)
% priv_inputSolution_Krylov - computes the set of input solutions in the Krylov
%    subspace
%
% Syntax:
%    linsys = priv_inputSolution_Krylov(linsys,options)
%
% Inputs:
%    linsys - linearSys object
%
% Outputs:
%    linsys - linearSys object
%    options - options for the computation of the reachable set
%
% Example:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Maximilian Perschl
% Written:       19-December-2016
% Last update:   25-April-2025 (MP, major refactor)
% Last revision: 21-March-2025 (TL)

% ------------------------------ BEGIN CODE -------------------------------

U = params.U;

% initialize variables
U_G = sparse(U.G);
nrOfGens = size(U_G,2);

% equivalent initial state
eqivState = [sparse([],[],[],linsys.nrOfDims,1);1];

% init Krylov order
KrylovOrder = 1;


% init
Ugen = [];
errorTaylor_sum_rad = sparse([],[],[],linsys.nrOfDims,1);

% consider generators; first check if generator is zero
for iGen = 1:nrOfGens
    U_g = U_G(:,iGen);
    if ~all(withinTol(U_g,0))
        
        % create integrator system; minus sign to fit paper
        A_int_g = [linsys.A, U_g; zeros(1,linsys.nrOfDims), 0];

        A_int_g = sparse(A_int_g);

        % Arnoldi
        [V_Ug,H_Ug,KrylovOrder] = ...
            priv_subspace_Krylov_individual_Jawecki(A_int_g,eqivState,KrylovOrder,options);

        % generate reduced order system for the generator
        A_red_g = H_Ug;
        B_red_g = V_Ug;
        %initialize linear reduced dynamics
        linRedSys_g = linearSys('linearReducedDynamics',A_red_g,B_red_g');
        linRedSys_g.taylor = taylorLinSys(A_red_g);
        
        % Loop through Taylor terms
        linRedSys_g.taylor = taylorLinSys(A_red_g);
        for i=1:options.taylorTerms+1
            % unprojected, partial result
            Apower_i = getTaylor(linRedSys_g,'Apower',struct('ithpower',i));
            U_unprojected = V_Ug*Apower_i(:,1)*options.timeStep^i/factorial(i);
            % compute sums
            Ugen(:,end+1) = U_unprojected(1:linsys.nrOfDims); 
        end

        % compute exponential matrix
        E = priv_expmRemainder(linRedSys_g,options.timeStep,options.taylorTerms);
        
        % error due to finite Taylor series
        errorTaylor_g = supremum(V_Ug*E(:,1)); 
        errorTaylor_sum_rad = errorTaylor_sum_rad + errorTaylor_g(1:linsys.nrOfDims);
    end
    % iGen
end

% round to 0 for numerical stability
errorTaylor_sum_rad(abs(errorTaylor_sum_rad) < eps) = 0;
% sparsify
errorTaylor_sum_rad = sparse(errorTaylor_sum_rad);
errorTaylor_sum = interval(-errorTaylor_sum_rad,errorTaylor_sum_rad);

% input zonotope without error terms
inputSolV = zonotope(sparse([],[],[],linsys.nrOfDims,1),Ugen);

% Check if error due to finite Taylor series needs to be considered
if options.krylovError > 2*eps
    % error due to order reduction
    errorRed = interval(-ones(linsys.nrOfDims,1),ones(linsys.nrOfDims,1)) * ...
        options.krylovError*size(U.G,2);
    % total error
    error_set = errorTaylor_sum + errorRed;
else
    error_set = errorTaylor_sum;
end

% initialize error for adaptive algorithm
linsys.krylov.total_U_0_error = radius(error_set);

% final input set
inputSolV = inputSolV + zonotope(error_set);


%write to object structure
linsys.krylov.V = U;
linsys.krylov.RV = inputSolV;

C = linsys.C;

% if clause in case this method is called again in adaptive case
if ~isfield(linsys.krylov,'Rpar_proj')
    linsys.krylov.Rpar_proj = zeros(linsys.nrOfOutputs,1);%interval(C*inputSolV);
end
linsys.krylov.Rpar_proj_0 = C*inputSolV;

linsys.krylov.RV_0 = inputSolV;

end

% ------------------------------ END OF CODE ------------------------------
