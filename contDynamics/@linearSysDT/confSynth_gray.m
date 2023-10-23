function [params, union_y_a, obj] = confSynth_gray(obj,params,options)
% confSynth_gray - specific conformance synthesis method for linear
%   discrete-time systems. This method identifies parameters of a gray-box 
%   model using fmincon in an outer loop as in Sec. III.B of [1] and 
%   updates the required uncertain sets according to [2],[3] in an inner 
%   loop.
%
% Syntax:
%    [params, union_y_a] = confSynth_gray(obj,params,options)
%
% Inputs:
%    obj - discrete-time linear system object
%    params - parameters defining the conformance problem
%    options - options for the conformance checking
%
% Outputs:
%    params - parameters solving the conformance problem
%    union_y_a - unified test cases (only deviation to nominal solution)
%    obj - updated discrete-time linear system object
%
% References:
%    [1] S. B. Liu and M. Althoff, "Online Verification of 
%        Impact-Force-Limiting Control for Physical Human-Robot 
%        Interaction," 2021 IEEE/RSJ International Conference on 
%        Intelligent Robots and Systems (IROS), 2021, pp. 777-783.
%    [2] M. Althoff. Checking and Establishing Reachset Conformance in
%        CORA 2023. In Proc. of the 10th International Workshop on 
%        Applied Verification of Continuous and Hybrid Systems, 2023.
%    [3] Liu et al., "Guarantees for Real Robotic Systems: Unifying Formal
%        Controller Synthesis and Reachset-Conformant Identification", 202x.
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Stefan Liu
% Written:       29-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Get variables from the parameter vector
param0 = options.p_0;
free   = options.p_free;
lb     = options.p_min;
ub     = options.p_max;

% Check if lb and ub exists
if any(isinf(lb)) || any(isinf(ub))
    error('Upper and lower bounds MUST be specified for variables.')
end

% Set options for fmincon
opts = optimoptions('fmincon','MaxFunEvals',10000,'MaxIter',10000,'Display','iter','Algorithm','sqp','UseParallel',false);

% Normalize variables for optimization and remove non-tunable params
% Why: https://stackoverflow.com/questions/12356041/is-normalization-useful-necessary-in-optimization
x0  = param0(free);
x0_ = (x0-lb(free))./(ub(free)-lb(free));
lb_ = zeros(size(x0));
ub_ = ones(size(x0));

% Declare variables for optimal params and corresponding cost
xLast     = []; % Last place computeall was called

% Call fmincon
[param_opt_, ~, ~, ~] = fmincon(@nest_objfun, x0_, [], [], [], [], lb_, ub_, [], opts); %nest_objfun nested below

% Nested objective function 
function cost = nest_objfun(x_)
    if ~isequal(x_,xLast) % Check if computation is necessary

        % Inner identification works on unormalized data 
        % -> hence unnormalize
        x = x_ .* (ub(free) - lb(free)) + lb(free);

        % Combine tunable and non-tunable parameters
        param0(free) = x;

        % Get the system matrices from the parameter vector
        [obj.A, obj.B, obj.C, obj.D] = params.getMatricesFromP(param0,obj);

        % Try the white identification method. This sometimes result in
        % an error when the new estimated matrices make the estimation
        % infeasible. We catch these errors and set the corresponding
        % cost to inf. 
        try
            [~, ~, cost] = obj.confSynth_dyn(params, options); % Just interested in the cost
        catch err
            % Display the error message with additional information
            disp('Error occurred during the identification:');
            disp(['Error message: ', err.message]);

            % Given matrices are infeasible, set cost to inf
            cost = inf;
        end

        xLast = x_;
    end
end


% Unnormalize data
param_opt = param_opt_ .* (ub(free) - lb(free)) + lb(free);

% Get the system matrices from the parameter vector 
param0(free) = param_opt;   
[obj.A, obj.B, obj.C, obj.D] = params.getMatricesFromP(param0,obj);

% Perform one last synthesis to obtain the optimal parameters
[params, union_y_a] = obj.confSynth_dyn(params, options);

end

% ------------------------------ END OF CODE ------------------------------
