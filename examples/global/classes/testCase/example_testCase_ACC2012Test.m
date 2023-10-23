function res = example_testCase_ACC2012Test()
% example_testCase_ACC2012Test - Creates test cases (see Def. 7 in [1])
% for the conformance checking procedure in [2].
%
% Syntax:
%    obj = example_testCase_ACC2012Test()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% References:
%    [1] Liu et al., "Guarantees for Real Robotic Systems: Unifying Formal
%        Controller Synthesis and Reachset-Conformant Identification", 2022
%    [2] M. Althoff and J. M. Dolan. Reachability computation of low-order 
%        models for the safety verification of high-order road vehicle 
%        models. In Proc. of the American Control Conference, 
%        page 3559–3566, 2012.
%
% Example:
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Matthias Althoff
% Written:       14-June-2023             
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% initialize result
ACC2012Test = cell(3,1);

%% loop over test cases
for iCase = 1:3
    % select test maneuver
    switch iCase
        case 1 % evasive test
            trajOptions = aux_man_evasiveTest(); 
            name = 'evasive test';
        case 2 % moose test
            trajOptions = aux_man_mooseTest();
            name = 'moose test';
        case 3 % cornering test
            trajOptions = aux_man_cornering();
            name = 'cornering test';
    end
    
    %% set options 
    params.tStart = 0; %start time
    params.x0 = [0; 0; 0; 15; 0; 0; 0]; %initial state for simulation
    params.timeStep = 0.01; %time step size for reachable set computation
 
    %reference trajectory
    [v_d, Psi_d, Psi_dot_d, x_d, y_d, ~, finalTime] = aux_vehicleTrajectory_smooth(trajOptions, params.timeStep, params.x0);

    % plot(x_d,y_d);
    % figure
    % plot(v_d);
    
    % remove values from previous iteration
    u = [];

    % create input trajectory
    u(1,:) = x_d; %x-position
    u(2,:) = y_d; %y-position
    u(3,:) = Psi_d; %orientation
    u(4,:) = Psi_dot_d; %yaw rate
    u(5,:) = v_d; %velocity

    % no sensor noise considered
    y = 0*u;
    
    % no disturbance added
    dim = 7; % system dimension
    w = zeros(7,size(u,2));

    % reference vector contains desired values and sensor noise
    uVec = [u;y;w];

    %specify continuous dynamics
    model = @vehicleDynamics_ST_controlled_BMW;
    carDyn = nonlinearSys(model,dim,size(uVec,1));

    %generate ode options
    stepsizeOptions = odeset('MaxStep',0.1*params.timeStep);
    options = odeset(stepsizeOptions);

    %compute single simulation
    inputChanges = ceil(finalTime/params.timeStep);
    xStep(:,1) = params.x0;
    params.tFinal = 0; % final time changed to start partial simulations

    for iChange = 1:inputChanges
        %reset options
        params.tStart = params.tFinal;
        params.tFinal = params.tFinal+finalTime/inputChanges;
        params.u = uVec(:,iChange);
        if iChange > 1
            params.x0 = x(end,:);
        end

        %simulate hybrid automaton
        [~,x] = simulate(carDyn, params, options); 
        xStep(:,iChange+1) = x(end,:); 
    end

    % save in test suite
    ACC2012Test{iCase} = testCase(xStep', uVec', xStep', params.timeStep, name, model); % measurement equals state in these test cases
end

%%save test suite
% path
path = [CORAROOT filesep 'models' filesep 'testCases' filesep 'autonomousDriving'];
if ~isfolder(path)
    mkdir([CORAROOT filesep 'models' filesep 'testCases'],'autonomousDriving');
end
% save
save([path filesep 'ACC2012Test'], 'ACC2012Test');

% example completed
res = true;

end


% Auxiliary functions -----------------------------------------------------

%% creation of smooth trajectory for desored values
function [v_d, Psi_d, Psi_dot_d, x_d, y_d, aVec_x, t_final] = aux_vehicleTrajectory_smooth(options,r,x0)

    %get data from options struct
    j_max = options.j_max; %[m/s^3]
    t = options.t;
    a_goal = options.a_goal;
    a_dir = options.a_dir;

    %obtain longitudinal and lateral acceleration
    for i=1:length(a_goal)
        a_goal_x(i) = cos(a_dir(i))*a_goal(i);
        a_goal_y(i) = sin(a_dir(i))*a_goal(i);
    end

    %initialize aVec_x and aVec_y
    aVec_x = a_goal_x(1);
    aVec_y = a_goal_y(1);

    for iPhase = 1:length(t)

        %determine difference in acceleration
        if iPhase == 1
            delta_a_x = 0;
            delta_a_y = 0;
        else
            delta_a_x = a_goal_x(iPhase) - a_goal_x(iPhase-1);
            delta_a_y = a_goal_y(iPhase) - a_goal_y(iPhase-1);
        end

        %normalize
        delta_abs = sqrt(delta_a_x^2 + delta_a_y^2);
        if delta_abs == 0
            delta_a_x = 0;
            delta_a_y = 0; 
        else
            delta_a_x = j_max*r*delta_a_x/delta_abs;
            delta_a_y = j_max*r*delta_a_y/delta_abs;
        end

        for iStep = 1:ceil(t(iPhase)/r)
            %new accelerations
            a_x_new = aVec_x(end) + delta_a_x;
            a_y_new = aVec_y(end) + delta_a_y;

            if abs(delta_a_x) > abs(delta_a_y)
                a_extr = a_goal_x(iPhase);
                a_comp = a_x_new;
                a_sign = sign(delta_a_x);
            else
                a_extr = a_goal_y(iPhase);
                a_comp = a_y_new;
                a_sign = sign(delta_a_y);
            end

            if (a_comp > a_extr) && (a_sign > 0)
                a_x_new = a_goal_x(iPhase);
                a_y_new = a_goal_y(iPhase);
            elseif (a_comp < a_extr) && (a_sign < 0)
                a_x_new = a_goal_x(iPhase);
                a_y_new = a_goal_y(iPhase);
            end

            %include new acceleration
            aVec_x(end+1) = a_x_new;
            aVec_y(end+1) = a_y_new;
        end
    end

    %auxiliary extension
    aVec_x(end+1) = aVec_x(end);
    aVec_y(end+1) = aVec_y(end);


    %time steps
    timeSteps = length(aVec_x)-1;

    %initialize
    Psi_d(1) = x0(2); %[rad]
    Psi_dot_d(1) = x0(3); %[rad/s]
    v_d(1) = x0(4); %[m/s]
    x_d(1) = x0(5);
    y_d(1) = x0(6);

    %compute Psi_d, y_d, a_x_d
    for iStep = 1:timeSteps
        %update velocity
        v_d(iStep+1) = v_d(iStep) + aVec_x(iStep)*r;
        %update Psi_d
        Psi_dot_d(iStep+1) = aVec_y(iStep+1)/v_d(iStep+1);
        Psi_d(iStep+1) = Psi_d(iStep) + Psi_dot_d(iStep)*r;
        %update x_d
        x_dot = cos(Psi_d(iStep))*v_d(iStep);
        x_d(iStep+1) = x_d(iStep) + x_dot*r;
        %y_d
        y_dot = sin(Psi_d(iStep))*v_d(iStep);
        y_d(iStep+1) = y_d(iStep) + y_dot*r;
    end

    %get t_final
    t_final = timeSteps*r;

end


%% encoding of cornering test case
function options = aux_man_cornering()

    %define maximum acceleration
    a_max = 6; %[m/s^2]

    %define maximum jerk
    j_max = 50; %[m/s^3]

    %set lateral and longitudinal acceleration
    t(1) = 0.4;
    t(2) = 1;
    t(3) = 1;
    t(4) = 0.4;


    %define a_goal
    a_goal(1) = 0;
    a_goal(2) = a_max;
    a_goal(3) = 0.8*a_max;
    a_goal(4) = 0;


    %define dir
    a_dir(1) = 0;
    a_dir(2) = +0.7*pi;
    a_dir(3) = +0.3*pi;
    a_dir(4) = 0;


    %write in options structure
    options.t = t;
    options.a_goal = a_goal;
    options.a_dir = a_dir;
    options.j_max = j_max;

end


%% encoding of most test case
function options = aux_man_mooseTest()

    %define maximum acceleration
    a_max = 6; %[m/s^2]

    %define maximum jerk
    j_max = 50; %[m/s^3]

    %set lateral and longitudinal acceleration
    t(1) = 0.4;
    t(2) = 1 - a_max/j_max;
    t(3) = 1;
    t(4) = 1;
    t(5) = 1 - a_max/j_max;
    t(6) = 1;
    t(7) = 0.4;

    %define a_goal
    a_goal(1) = 0;
    a_goal(2) = a_max;
    a_goal(3) = a_max;
    a_goal(4) = 0;
    a_goal(5) = a_max;
    a_goal(6) = a_max;
    a_goal(7) = 0;

    %define dir
    a_dir(1) = 0;
    a_dir(2) = +0.5*pi;
    a_dir(3) = -0.5*pi;
    a_dir(4) = 0;
    a_dir(5) = -0.5*pi;
    a_dir(6) = +0.5*pi;
    a_dir(7) = 0;

    %write in options structure
    options.t = t;
    options.a_goal = a_goal;
    options.a_dir = a_dir;
    options.j_max = j_max;

end

%% encoding of evasive test case
function options = aux_man_evasiveTest()

    %define maximum acceleration
    a_max = 6; %[m/s^2]

    %define maximum jerk
    j_max = 50; %[m/s^3]

    %set lateral and longitudinal acceleration
    t(1) = 0.4;
    t(2) = 0.75;
    t(3) = 0.63;
    t(4) = 0.65;


    %define a_goal
    a_goal(1) = 0;
    a_goal(2) = a_max;
    a_goal(3) = a_max;
    a_goal(4) = 0;


    %define dir
    a_dir(1) = 0;
    a_dir(2) = +0.75*pi;
    a_dir(3) = -0.75*pi;
    a_dir(4) = -pi;


    %write in options structure
    options.t = t;
    options.a_goal = a_goal;
    options.a_dir = a_dir;
    options.j_max = j_max;

end

% ------------------------------ END OF CODE ------------------------------
