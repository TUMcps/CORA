function R = priv_observe_intersectionFreeAdaptive(nlnsysDT,params,options)
% priv_observe_intersectionFreeAdaptive - computes the guaranteed state
%    estimation approach according to the intersection-free approach
%    when the gain changes in each iteration, see [1], [2]; the approach is 
%    extended here for nonlinear systems.
%
% Syntax:
%    R = priv_observe_intersectionFreeAdaptive(nlnsysDT,params,options)
%
% Inputs:
%    nlnsysDT - nonlinearSysDT object
%    params - model parameters
%    options - options for the guaranteed state estimation
%
% Outputs:
%    R - observed set of points in time
%
% Reference:
%    [1] Christophe Combastel. Zonotopes and kalman observers:
%        Gain optimality under distinct uncertainty paradigms and
%        robust convergence. Automatica, 55:265-273, 2015.
%    [2] M. Althoff and J. J. Rath. Comparison of Set-Based Techniques 
%        for Guaranteed State Estimation of Linear Disturbed Systems, 
%        in preparation.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       25-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%time period
tVec = params.tStart:options.timeStep:params.tFinal-options.timeStep;

% initialize parameter for the output equation
R = cell(length(tVec),1);
% store first reachable set
Rnext.tp = params.R0;
R{1} = params.R0;

% F is chosen in [1] such that they are multiplied with unit
% uncertainties; thus, E and F can be seen as generators of zonotopes
% representing the disturbance and noise set
F = generators(params.V);

% loop over all time steps
for k = 1:length(tVec)-1
    
    %% linearize nonlinear system
    %linearization point p
    p.u = params.uTransVec(:,k);
    p.x = center(R{k});

    %substitute p into the system equation in order to obtain the constant
    %input
    f0 = nlnsysDT.mFile(p.x, p.u);

    %get jacobian matrices
    [A_lin,B_lin] = nlnsysDT.jacobian(p.x, p.u);

    %save linearization point
    nlnsysDT.linError.p=p;

    %translate Rinit by linearization point
    Rdelta = R{k} + (-nlnsysDT.linError.p.x);

    % obtain linearization error
    if options.tensorOrder == 2
        Verror = priv_linError_mixed_noInt(nlnsysDT, Rdelta, params, options);
    elseif options.tensorOrder == 3
        Verror = priv_linError_thirdOrder(nlnsysDT, Rdelta, params, options);
    end
    
    % Compute observer gain
    if options.observerType == 1 % Combastel, FRad-C in [2]
        % obtain generators
        G = generators(Rnext.tp);
        G_comb = G*G';
        L = A_lin*G_comb*nlnsysDT.C'*inv(nlnsysDT.C*G_comb*nlnsysDT.C' + F*F');
    else
        disp('this observer type is not yet implemented')
    end
    
    % Prediction, eq. (11) in [2]
    Rnext.tp = (A_lin-L*nlnsysDT.C)*Rdelta - L*nlnsysDT.C*nlnsysDT.linError.p.x  + ...
        L*params.y(:,k) + (-L)*params.V + params.W + Verror + f0;
    
    % Order reduction
    Rnext.tp = reduce(Rnext.tp,options.reductionTechnique,options.zonotopeOrder);
 
    % Store result
    R{k+1} = Rnext.tp;
end

% ------------------------------ END OF CODE ------------------------------
