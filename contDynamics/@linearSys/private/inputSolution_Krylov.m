function [obj] = inputSolution_Krylov(obj, options)
% inputSolution_Krylov - computes the set of input solutions in the Krylov
% subspace
%
% Syntax:  
%    [obj] = inputSolution_Krylov(obj,options)
%
% Inputs:
%    obj - linearSys object
%
% Outputs:
%    obj - linearSys object
%    options - options for the computation of the reachable set
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      19-December-2016
% Last update:  28-October-2017
%               25-October-2018
% Last revision:---

%------------- BEGIN CODE --------------


% set of possible inputs
V = obj.B*options.U;

% separate input set into center and generators
U_c = sparse(center(V));
U_G = sparse(generators(V));
nrOfGens = length(U_G(1,:));

% equivalent initial state
eqivState = [zeros(length(obj.A),1); 1];

% step size
delta = options.timeStep;

% init Krylov order
KrylovOrder = 1;


% concider center; first check if center is zero
if ~all(U_c == 0)
    % obtain constant input solution
    [Usum, errorTaylor_c, ~, KrylovOrderConst] = ...
        constInputSolution(obj, U_c, KrylovOrder, options, 0);
else
    Usum = zeros(length(obj.A),1);
    errorTaylor_c = zeros(length(obj.A),1);
    KrylovOrderConst = 1;
end

% init
Ugen = [];
errorTaylor_sum = errorTaylor_c;

% concider generators; first check if generator is zero
for iGen = 1:nrOfGens
    U_g = U_G(:,iGen);
    if ~all(U_g == 0)
        
        % create integrator system; minus sign to fit paper
        A_int_g = [obj.A, U_g; zeros(1,length(obj.A)), 0];

        % minimum eigenvalue
        %nu_A_int_g = eigs((A_int_g+A_int_g')/2, 1, 'sa'); % lambda_min(A+A*/2) %<-- for Wang approx.
        nu_A_int_g = eigs((A_int_g+A_int_g')/2, 1, 'lm'); % lambda_max(A+A*/2) %<-- for Jia approx.

        % Arnoldi
        [V_Ug,H_Ug,KrylovOrder] = subspace_Krylov_individual(A_int_g,nu_A_int_g,eqivState,KrylovOrder,options);

        % change results back \dot{x} = Ax (Krylov papers use \dot{x} = -Ax);
        % V remains unchanged
        %H_Ug = -H_Ug;
        
        % generate reduced order system for the generators
        A_red_g = H_Ug;
        B_red_g = V_Ug;
        %initialize linear reduced dynamics
        linRedSys_g = linearSys('linearReducedDynamics',A_red_g,B_red_g);

        % compute exponential matrix
        linRedSys_g = exponential(linRedSys_g,options);

        % Loop through Taylor terms
        for i=1:options.taylorTerms
            % unprojected, partial result
            U_unprojected = V_Ug*linRedSys_g.taylor.powers{i}(:,1)*delta^i/factorial(i);
            % compute sums
            Ugen(:,end+1) = U_unprojected(1:length(obj.A)); 
        end
        
        % error due to finite Taylor series
        errorTaylor_g = V_Ug*linRedSys_g.taylor.error(:,1); 
        errorTaylor_sum = errorTaylor_sum + errorTaylor_g(1:length(obj.A));
    end
    iGen
end

%compute vTrans 
if iscell(options.uTrans)
    timeSteps = length(options.uTrans(1,:)); % number of time steps
    vTrans = cell(timeSteps); % init vTrans
    for i=1:timeSteps
        vTrans{i}=obj.B*options.uTrans{i};
    end
else
    vTrans=obj.B*options.uTrans;
end

if ~all(vTrans == 0) % input trajectory not yet implemented
    
    %compute additional uncertainty if origin is not contained in input set
    if options.originContained
        [inputSolVtrans, errorTaylor_vTrans] = ...
            constInputSolution(obj, vTrans, KrylovOrderConst, options, 0);
        inputCorr=zeros(dim,1);
    else
        % obtain constant input solution
        [inputSolVtrans, errorTaylor_vTrans, inputCorr] = ...
            constInputSolution(obj, vTrans, KrylovOrderConst, options, 1);
    end
else
    inputSolVtrans = zeros(length(obj.A),1);
    inputCorr = zeros(length(obj.A),1);
    errorTaylor_vTrans = zeros(length(obj.A),1);
end

% input zonotope without error terms
inputSolV = zonotope([Usum,Ugen]);

% Check if error due to finite Taylor series needs to be considered
if options.krylovError > 2*eps
    % error due to order reduction
    errorRed = interval(-ones(length(obj.A),1),ones(length(obj.A),1)) * ...
        options.krylovError*norm(V+vTrans)*delta;
    % total error
    error = errorTaylor_sum + errorRed + errorTaylor_vTrans;
    % final input set
    inputSolV = inputSolV + zonotope(error);
end


%write to object structure
obj.taylor.V = V;
obj.taylor.RV = inputSolV;
obj.taylor.Rtrans = inputSolVtrans; % to avoid additional generators, all errors are added to inputSolV
obj.taylor.inputCorr = inputCorr;
obj.taylor.eAtInt = [];

end


function [Usum, errorTaylor, inputCorr, KrylovOrder] = ...
    constInputSolution(obj, input, KrylovOrder, options, inputTieFlag)

    inputCorr = [];

    % equivalent initial state
    eqivState = [zeros(length(obj.A),1); 1];
    
    % init Usum
    Usum = zeros(length(obj.A),1);

    % create integrator system
    A_int = [obj.A, input; zeros(1,length(obj.A)), 0];

    % minimum eigenvalue
    %nu_A_int_g = eigs((A_int_g+A_int_g')/2, 1, 'sa'); % lambda_min(A+A*/2) %<-- for Wang approx.
    nu_A_int = eigs((A_int+A_int')/2, 1, 'lm'); % lambda_max(A+A*/2) %<-- for Jia approx.
    
    % Arnoldi
    [V_U,H_U,KrylovOrder] = subspace_Krylov_individual(A_int,nu_A_int,eqivState,KrylovOrder,options);
    
    % change results back \dot{x} = Ax (Krylov papers use \dot{x} = -Ax);
    % V remains unchanged
    %H_U = -H_U;
    
    % generate reduced order system for the generators
    A_red = H_U;
    B_red = V_U;
    %initialize linear reduced dynamics
    linRedSys = linearSys('linearReducedDynamics',A_red,B_red);

    % compute exponential matrix
    linRedSys = exponential(linRedSys, options);
    
    % Loop through Taylor terms
    for i=1:options.taylorTerms
        % unprojected, partial result
        U_unprojected = V_U*linRedSys.taylor.powers{i}(:,1) * ...
            options.timeStep^i/factorial(i);
        % compute sums
        Usum = Usum + U_unprojected(1:length(obj.A)); 
    end
    
    % error due to finite Taylor series
    errorTaylor_unproj = V_U*linRedSys.taylor.error(:,1); 
    errorTaylor = errorTaylor_unproj(1:length(obj.A));
    
    % if inputTie is required
    if inputTieFlag
        % compute time interval error (tie) for inputs
        linRedSys = inputTie(linRedSys,options);
        F_tmp = linRedSys.taylor.inputF;
        %tie_error 
        Fmid = center(F_tmp);
        Frad = rad(F_tmp);
        inputTie_error_mid = V_U*Fmid(:,1); % norm(equiv_state) = 1;
        inputTie_error_rad = abs(V_U*Frad(:,1)); % norm(equiv_state) = 1;
        inputCorr_unprojected = interval(inputTie_error_mid - ...
            inputTie_error_rad, inputTie_error_mid + inputTie_error_rad); 
        inputCorr = zonotope(inputCorr_unprojected(1:length(obj.A))); 
    end
end


%------------- END OF CODE --------------