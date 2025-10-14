function res = test_hybridAutomaton_identify
% test_hybridAutomaton_identify - unit test for system identification
%
% Syntax:
%    res = test_hybridAutomaton_identify
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Niklas Kochdumper
% Written:       13-June-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % Test 1: system without inputs, but with jumps
    HAorig = bouncing_ball(-0.75);

    simOpts.x0 = [2;0];
    simOpts.tFinal = 3;
    simOpts.startLoc = 1;

    % simulate system
    [t_1,x_1] = simulate(HAorig,simOpts);

    simOpts.x0 = [1.5;0.1];
    simOpts.tFinal = 2.5;
    
    % simulate system
    [t_2,x_2] = simulate(HAorig,simOpts);

    HA = hybridAutomaton.identify(x_1,t_1);
    assert(isa(HA,'hybridAutomaton'));

    t{1} = t_1; t{2} = t_2;
    x{1} = x_1; x{2} = x_2;
    HA = hybridAutomaton.identify(x,t);
    assert(isa(HA,'hybridAutomaton'));

    % identify using trajectory object
    traj_1 = trajectory([],x_1,[],t_1);
    HA = hybridAutomaton.identify(traj_1);
    assert(isa(HA,'hybridAutomaton'));

    traj = [traj_1; trajectory([],x_2,[],t_2)];
    HA = hybridAutomaton.identify(traj);
    assert(isa(HA,'hybridAutomaton'));

    % Test 2: system with inputs
    HAorig = roomHeating();

    simOpts.x0 = [5;19];
    simOpts.tFinal = 5;
    simOpts.u = 20;
    simOpts.startLoc = 1;

    % simulate system
    [t_1,x_1] = simulate(HAorig,simOpts);
    u_1 = simOpts.u*ones(1,length(t_1));

    simOpts.x0 = [18,10];
    simOpts.tFinal = 6;
    simOpts.u = 19;

    % simulate system
    [t_2,x_2] = simulate(HAorig,simOpts);
    u_2 = simOpts.u*ones(1,length(t_2));

    HA = hybridAutomaton.identify(x_1,t_1,u_1);
    assert(isa(HA,'hybridAutomaton'));

    t{1} = t_1; t{2} = t_2;
    x{1} = x_1; x{2} = x_2;
    u{1} = u_1; u{2} = u_2;
    HA = hybridAutomaton.identify(x,t,u);
    assert(isa(HA,'hybridAutomaton'));

    % identify using trajectory object
    traj_1 = trajectory(u_1,x_1,[],t_1);
    HA = hybridAutomaton.identify(traj_1);
    assert(isa(HA,'hybridAutomaton'));

    traj = [traj_1; trajectory(u_2,x_2,[],t_2)];
    HA = hybridAutomaton.identify(traj);
    assert(isa(HA,'hybridAutomaton'));

    % test successfull if it runs through without errors
    res = true;

% ------------------------------ END OF CODE ------------------------------
