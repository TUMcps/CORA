function R = priv_observe_interval(linsysDT,params,options)
% priv_observe_interval - computes the guaranteed state estimation 
%    approach according to the interval approach, see Alg. 1 in [1].
%
% Syntax:
%    R = priv_observe_interval(linsysDT,params,options)
%
% Inputs:
%    linsysDT - discrete-time linear system object
%    params - model parameters
%    options - options for the guaranteed state estimation
%
% Outputs:
%    R - observed set of points in time
%
% Reference:
%    [1] M. Althoff and J. J. Rath. Comparison of Set-Based Techniques 
%        for Guaranteed State Estimation of Linear Disturbed Systems, 
%        in preparation.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       05-January-2021
% Last update:   25-February-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%%initialize computation

%time period
tVec = params.tStart:options.timeStep:params.tFinal-options.timeStep;
timeSteps = length(tVec);

% initialize parameter for the output equation
R = cell(length(tVec),1);

% store first reachable set
R{1} = params.R0;

%% initialize
x = center(params.R0);
D = -options.L*params.V + params.W;
R_x = params.R0;
R_wv = zeros(length(x),1);

%% loop over all time steps
for i = 1:timeSteps-1
    
    % error set
    E = box(R_x) + R_wv;
    
    % reachable set
    R{i+1} = x + E;
    
    % update state
    x = linsysDT.A*x + linsysDT.B*params.uTransVec(:,i) ...
        + options.L*(params.y(:,i) - linsysDT.C*x);

    % Update auxiliary sets
    R_x = (linsysDT.A - options.L*linsysDT.C)*R_x; 
    R_wv = R_wv + box(D);
    D = (linsysDT.A - options.L*linsysDT.C)*D;
end

% ------------------------------ END OF CODE ------------------------------
