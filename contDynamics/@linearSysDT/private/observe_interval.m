function R = observe_interval(obj,options)
% observe_interval - computes the guaranteed state estimation 
% approach according to the interval approach, see Alg. 1 in [1].
%
%
% Syntax:
%    R = observe_interval(obj,options)
%
% Inputs:
%    obj - discrete-time linear system object
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
% Example: 
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
tVec = options.tStart:options.timeStep:options.tFinal-options.timeStep;
timeSteps = length(tVec);

% initialize parameter for the output equation
R = cell(length(tVec),1);

% store first reachable set
R{1} = options.R0;


%% initialize
x = center(options.R0);
D = -options.L*options.V + options.W;
R_x = options.R0;
R_wv = zeros(length(x),1);


%% loop over all time steps
for i = 1:timeSteps-1
    
    % error set
    E = box(R_x) + R_wv;
    
    % reachable set
    R{i+1} = x + E;
    
    % update state
    x = obj.A*x + obj.B*options.uTransVec(:,i) + options.L*(options.y(:,i) - obj.C*x);

    % Update auxiliary sets
    R_x = (obj.A - options.L*obj.C)*R_x; 
    R_wv = R_wv + box(D);
    D = (obj.A - options.L*obj.C)*D;
end


% ------------------------------ END OF CODE ------------------------------
