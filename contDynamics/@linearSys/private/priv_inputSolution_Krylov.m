function linsys = priv_inputSolution_Krylov(linsys,options)
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

% Authors:       Matthias Althoff
% Written:       19-December-2016
% Last update:   28-October-2017
%                25-October-2018
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set of possible inputs
V = linsys.B*options.U;

% separate input set into center and generators
U_c = sparse(V.c);
U_G = sparse(V.G);
nrOfGens = size(U_G,2);

% equivalent initial state
eqivState = [zeros(linsys.nrOfDims,1); 1];

% init Krylov order
KrylovOrder = 1;

% consider center; first check if center is zero
if ~all(withinTol(U_c,0))
    % obtain constant input solution
    [Usum, errorTaylor_c, ~, KrylovOrderConst] = ...
        aux_constInputSolution(linsys, U_c, KrylovOrder, options, false);
else
    Usum = zeros(linsys.nrOfDims,1);
    errorTaylor_c = zeros(linsys.nrOfDims,1);
    KrylovOrderConst = 1;
end

% init
Ugen = [];
errorTaylor_sum = errorTaylor_c;

% consider generators; first check if generator is zero
for iGen = 1:nrOfGens
    U_g = U_G(:,iGen);
    if ~all(withinTol(U_g,0))
        
        % create integrator system; minus sign to fit paper
        A_int_g = [linsys.A, U_g; zeros(1,linsys.nrOfDims), 0];

        % minimum eigenvalue
        %nu_A_int_g = eigs((A_int_g+A_int_g')/2, 1, 'sa'); % lambda_min(A+A*/2) %<-- for Wang approx.
        nu_A_int_g = eigs((A_int_g+A_int_g')/2, 1, 'lm'); % lambda_max(A+A*/2) %<-- for Jia approx.

        % Arnoldi
        [V_Ug,H_Ug,KrylovOrder] = ...
            priv_subspace_Krylov_individual(A_int_g,nu_A_int_g,eqivState,KrylovOrder,options);

        % change results back \dot{x} = Ax (Krylov papers use \dot{x} = -Ax);
        % V remains unchanged
        %H_Ug = -H_Ug;
        
        % generate reduced order system for the generators
        A_red_g = H_Ug;
        B_red_g = V_Ug;
        %initialize linear reduced dynamics
        linRedSys_g = linearSys('linearReducedDynamics',A_red_g,B_red_g);

        % Loop through Taylor terms
        for i=1:options.taylorTerms
            % unprojected, partial result
            Apower_i = getTaylor(linRedSys_g,'Apower',struct('ithpower',i));
            U_unprojected = V_Ug*Apower_i(:,1)*options.timeStep^i/factorial(i);
            % compute sums
            Ugen(:,end+1) = U_unprojected(1:linsys.nrOfDims); 
        end

        % compute exponential matrix
        % priv_expmRemainder?
        E = new_expmRemainder(linRedSys_g,options.timeStep,options.taylorTerms);
        
        % error due to finite Taylor series
        errorTaylor_g = V_Ug*E(:,1); 
        errorTaylor_sum = errorTaylor_sum + errorTaylor_g(1:linsys.nrOfDims);
    end
    iGen
end

%compute vTrans 
if iscell(options.uTrans)
    timeSteps = length(options.uTrans(1,:)); % number of time steps
    vTrans = cell(timeSteps); % init vTrans
    for i=1:timeSteps
        vTrans{i}=linsys.B*options.uTrans{i};
    end
else
    vTrans=linsys.B*options.uTrans;
end

if ~all(withinTol(vTrans,0)) % input trajectory not yet implemented
    
    %compute additional uncertainty if origin is not contained in input set
    if options.originContained
        [inputSolVtrans, errorTaylor_vTrans] = ...
            aux_constInputSolution(linsys, vTrans, KrylovOrderConst, options, false);
        inputCorr = zeros(linsys.nrOfDims,1);
    else
        % obtain constant input solution
        [inputSolVtrans, errorTaylor_vTrans, inputCorr] = ...
            aux_constInputSolution(linsys, vTrans, KrylovOrderConst, options, true);
    end
else
    inputSolVtrans = zeros(linsys.nrOfDims,1);
    inputCorr = zeros(linsys.nrOfDims,1);
    errorTaylor_vTrans = zeros(linsys.nrOfDims,1);
end

% input zonotope without error terms
inputSolV = zonotope(Usum,Ugen);

% Check if error due to finite Taylor series needs to be considered
if options.krylovError > 2*eps
    % error due to order reduction
    errorRed = interval(-ones(linsys.nrOfDims,1),ones(linsys.nrOfDims,1)) * ...
        options.krylovError*norm(V+vTrans)*options.timeStep;
    % total error
    err = errorTaylor_sum + errorRed + errorTaylor_vTrans;
    % final input set
    inputSolV = inputSolV + zonotope(err);
end


%write to object structure
linsys.taylor.V = V;
linsys.taylor.RV = inputSolV;
% to avoid additional generators, all errors of inputSolVtrans are added to inputSolV
linsys.taylor.Rtrans = inputSolVtrans;
linsys.taylor.inputCorr = inputCorr;
linsys.taylor.eAtInt = [];

end


% Auxiliary functions -----------------------------------------------------

function [Usum, errorTaylor, inputCorr, KrylovOrder] = ...
    aux_constInputSolution(sys, input, KrylovOrder, options, inputTieFlag)

    inputCorr = [];

    % equivalent initial state
    eqivState = [zeros(sys.nrOfDims,1); 1];
    
    % init Usum
    Usum = zeros(sys.nrOfDims,1);

    % create integrator system
    A_int = [sys.A, input; zeros(1,sys.nrOfDims), 0];

    % minimum eigenvalue
    %nu_A_int_g = eigs((A_int_g+A_int_g')/2, 1, 'sa'); % lambda_min(A+A*/2) %<-- for Wang approx.
    nu_A_int = eigs((A_int+A_int')/2, 1, 'lm'); % lambda_max(A+A*/2) %<-- for Jia approx.
    
    % Arnoldi
    [V_U,H_U,KrylovOrder] = ...
        priv_subspace_Krylov_individual(A_int,nu_A_int,eqivState,KrylovOrder,options);
    
    % change results back \dot{x} = Ax (Krylov papers use \dot{x} = -Ax);
    % V remains unchanged
    %H_U = -H_U;
    
    % generate reduced order system for the generators
    A_red = H_U;
    B_red = V_U;
    %initialize linear reduced dynamics
    linRedSys = linearSys('linearReducedDynamics',A_red,B_red);

    % compute exponential matrix
    % priv_expmRemainder?
    E = new_expmRemainder(linRedSys,options.timeStep,options.taylorTerms);
    
    % Loop through Taylor terms
    for i=1:options.taylorTerms
        % unprojected, partial result
        Apower_i = getTaylor(linRedSys_g,'Apower',struct('ithpower',i));
        U_unprojected = V_U*Apower_i(:,1)*options.timeStep^i/factorial(i);
        % compute sums
        Usum = Usum + U_unprojected(1:sys.nrOfDims); 
    end
    
    % error due to finite Taylor series
    errorTaylor_unproj = V_U*E(:,1);
    errorTaylor = errorTaylor_unproj(1:sys.nrOfDims);
    
    % if inputTie is required
    if inputTieFlag
        % priv_curvatureInput?
        G = new_curvatureInput(linRedSys,options.timeStep,options.taylorTerms);
        Gmid = center(G);
        Grad = rad(G);
        inputTie_error_mid = V_U*Gmid(:,1); % norm(equiv_state) = 1;
        inputTie_error_rad = abs(V_U*Grad(:,1)); % norm(equiv_state) = 1;
        inputCorr_unprojected = interval(...
            inputTie_error_mid - inputTie_error_rad, ...
            inputTie_error_mid + inputTie_error_rad); 
        inputCorr = zonotope(inputCorr_unprojected(1:sys.nrOfDims)); 
    end
end


% ------------------------------ END OF CODE ------------------------------
