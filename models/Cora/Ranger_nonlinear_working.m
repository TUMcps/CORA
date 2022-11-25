function HA = Ranger_nonlinear_working(~)


%% Generated on 13-May-2022

%-------------Automaton created from Component 'RoboWorld'-----------------

%% Interface Specification:
%   This section clarifies the meaning of state, input & output dimensions
%   by showing their mapping to SpaceEx variable names. 

% Component 1 (RoboWorld.RangerRoboChart):
%  state x := [ls; as; timer; turnTimer]
%  input u := [uDummy]

% Component 2 (RoboWorld.EnvironmentMapping.Move_OperationMapping):
%  state x := [arg_x; arg_y; timer]
%  input u := [yaw; ls; as]

% Component 3 (RoboWorld.EnvironmentMapping.Environment.RobotMovement):
%  state x := [pos_x; pos_y; vel_x; vel_y; acc_x; acc_y; yaw; yawVel; yawAcc; t]
%  input u := [arg_x; arg_y]

% Component 4 (RoboWorld.EnvironmentMapping.Environment.CollisionDetection):
%  state x := [pos_x; pos_y]
%  input u := [uDummy]

% Component 5 (RoboWorld.EnvironmentMapping.Environment.Obstacle_InputEventMapping):
%  state x := [pos_x; pos_y]
%  input u := [uDummy]

%------------------Component RoboWorld.RangerRoboChart---------------------

%-------------------------State Moving_setArgs-----------------------------

%% equation:
%   ls' == 0 &
%   as' == 0 &
%   timer' == 1 &
%   turnTimer' == 1
dynA = ...
[0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0];
dynB = ...
[0;0;0;0];
dync = ...
[0;0;1;1];
dynamics = linearSys(dynA, dynB, dync);

%% equation:
%   timer <= 0
A = ...
[0,0,1,0];
b = ...
[0];
polyOpt = struct('A', A, 'b', b);
inv = mptPolytope(polyOpt);

trans = {};
%% equation:
%   ls := lv &
%   as := 0
resetA = ...
[0,0,0,0;0,0,0,0;0,0,1,0;0,0,0,1];
resetc = ...
[1;0;0;0];
reset = struct('A', resetA, 'c', resetc,'hasInput',0);

%% equation:
%   
A = ...
zeros([0,0]);
b = ...
zeros([0,1]);
polyOpt = struct('A', A, 'b', b);
guard = mptPolytope(polyOpt);

trans{1} = transition(guard, reset, 2,[], 4);

loc{1} = location('S1', inv, trans, dynamics);



%---------------------------State Moving_call------------------------------

%% equation:
%   ls' == 0 &
%   as' == 0 &
%   timer' == 1 &
%   turnTimer' == 1
dynA = ...
[0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0];
dynB = ...
[0;0;0;0];
dync = ...
[0;0;1;1];
dynamics = linearSys(dynA, dynB, dync);

%% equation:
%   timer <= 0
A = ...
[0,0,1,0];
b = ...
[0];
polyOpt = struct('A', A, 'b', b);
inv = mptPolytope(polyOpt);

trans = {};
%% equation:
%   no reset equation given
resetA = ...
[1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1];
resetc = ...
[0;0;0;0];
reset = struct('A', resetA, 'c', resetc,'hasInput',0);

%% equation:
%   
A = ...
zeros([0,0]);
b = ...
zeros([0,1]);
polyOpt = struct('A', A, 'b', b);
guard = mptPolytope(polyOpt);

trans{1} = transition(guard, reset, 3,"moveCall", 4);

loc{2} = location('S2', inv, trans, dynamics);



%-----------------------------State Moving---------------------------------

%% equation:
%   ls' == 0 &
%   as' == 0 &
%   timer' == 1 &
%   turnTimer' == 1
dynA = ...
[0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0];
dynB = ...
[0;0;0;0];
dync = ...
[0;0;1;1];
dynamics = linearSys(dynA, dynB, dync);

%% equation:
%   
A = ...
zeros([0,0]);
b = ...
zeros([0,1]);
polyOpt = struct('A', A, 'b', b);
inv = mptPolytope(polyOpt);

trans = {};
%% equation:
%   timer := 0
resetA = ...
[1,0,0,0;0,1,0,0;0,0,0,0;0,0,0,1];
resetc = ...
[0;0;0;0];
reset = struct('A', resetA, 'c', resetc,'hasInput',0);

%% equation:
%   
A = ...
zeros([0,0]);
b = ...
zeros([0,1]);
polyOpt = struct('A', A, 'b', b);
guard = mptPolytope(polyOpt);

trans{1} = transition(guard, reset, 4,"obstacle", 4);

loc{3} = location('S3', inv, trans, dynamics);



%-------------------------State Turning_setArgs----------------------------

%% equation:
%   ls' == 0 &
%   as' == 0 &
%   timer' == 1 &
%   turnTimer' == 1
dynA = ...
[0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0];
dynB = ...
[0;0;0;0];
dync = ...
[0;0;1;1];
dynamics = linearSys(dynA, dynB, dync);

%% equation:
%   timer <= 0
A = ...
[0,0,1,0];
b = ...
[0];
polyOpt = struct('A', A, 'b', b);
inv = mptPolytope(polyOpt);

trans = {};
%% equation:
%   ls := 0 &
%   as := av &
%   timer := 0
resetA = ...
[0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,1];
resetc = ...
[0;1.5707999999999999740651901447563;0;0];
reset = struct('A', resetA, 'c', resetc,'hasInput',0);

%% equation:
%   timer >= 0
c = [0;0;1;0];
d = 0;

guard = conHyperplane(c,d);

trans{1} = transition(guard, reset, 5,[]);

loc{4} = location('S4', inv, trans, dynamics);



%--------------------------State Turning_call------------------------------

%% equation:
%   ls' == 0 &
%   as' == 0 &
%   timer' == 1 &
%   turnTimer' == 1
dynA = ...
[0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0];
dynB = ...
[0;0;0;0];
dync = ...
[0;0;1;1];
dynamics = linearSys(dynA, dynB, dync);

%% equation:
%   timer <= 0
A = ...
[0,0,1,0];
b = ...
[0];
polyOpt = struct('A', A, 'b', b);
inv = mptPolytope(polyOpt);

trans = {};
%% equation:
%   timer := 0
resetA = ...
[1,0,0,0;0,1,0,0;0,0,0,0;0,0,0,1];
resetc = ...
[0;0;0;0];
reset = struct('A', resetA, 'c', resetc,'hasInput',0);

%% equation:
%   timer >= 0
c = [0;0;1;0];
d = 0;

guard = conHyperplane(c,d);

trans{1} = transition(guard, reset, 6,"moveCall");

loc{5} = location('S5', inv, trans, dynamics);



%-----------------------------State Turning--------------------------------

%% equation:
%   turnTimer' == 1 &
%   ls' == 0 &
%   as' == 0 &
%   timer' == 1
dynA = ...
[0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0];
dynB = ...
[0;0;0;0];
dync = ...
[0;0;1;1];
dynamics = linearSys(dynA, dynB, dync);

%% equation:
%   timer <= turnTime
A = ...
[0,0,1,0];
b = ...
[2];
polyOpt = struct('A', A, 'b', b);
inv = mptPolytope(polyOpt);

trans = {};
%% equation:
%   timer := 0
resetA = ...
[1,0,0,0;0,1,0,0;0,0,0,0;0,0,0,1];
resetc = ...
[0;0;0;0];
reset = struct('A', resetA, 'c', resetc,'hasInput',0);

%% equation:
%   timer >= turnTime
c = [0;0;1;0];
d = 2;

guard = conHyperplane(c,d);

trans{1} = transition(guard, reset, 1,[]);

loc{6} = location('S6', inv, trans, dynamics);



comp{1} = hybridAutomaton(loc);

iBinds{1} = [0 1];
%-----Component RoboWorld.EnvironmentMapping.Move_OperationMapping---------

%---------------------------State WaitForCall------------------------------

%% equation:
%   arg_x' == 0 &
%   arg_y' == 0 &
%   timer' == 1
dynA = ...
[0,0,0;0,0,0;0,0,0];
dynB = ...
[0,0,0;0,0,0;0,0,0];
dync = ...
[0;0;1];
dynamics = linearSys(dynA, dynB, dync);

%% equation:
%   
A = ...
zeros([0,0]);
b = ...
zeros([0,1]);
polyOpt = struct('A', A, 'b', b);
inv = mptPolytope(polyOpt);

trans = {};
%% equation:
%   timer := 0
resetA = ...
[1,0,0;0,1,0;0,0,0];
resetc = ...
[0;0;0];
reset = struct('A', resetA, 'c', resetc,'hasInput',0);

%% equation:
%   
A = ...
zeros([0,0]);
b = ...
zeros([0,1]);
polyOpt = struct('A', A, 'b', b);
guard = mptPolytope(polyOpt);

trans{1} = transition(guard, reset, 2,"moveCall", 3);

loc{1} = location('S1', inv, trans, dynamics);



%-----------------------State SetLinearSpeedArgs---------------------------

%% equation:
%   arg_x' == 0 &
%   arg_y' == 0 &
%   timer' == 1
dynA = ...
[0,0,0;0,0,0;0,0,0];
dynB = ...
[0,0,0;0,0,0;0,0,0];
dync = ...
[0;0;1];
dynamics = linearSys(dynA, dynB, dync);

%% equation:
%   timer <= 0
A = ...
[0,0,1];
b = ...
[0];
polyOpt = struct('A', A, 'b', b);
inv = mptPolytope(polyOpt);

trans = {};
%% equation:
%   arg_x := ls*(-1 + 0.5*(yaw-3.14159)^2 - 0.041667*(yaw-3.14159)^4 + 0.0013889*(yaw-3.14159)^6 - 2.4801e-5*(yaw-3.14159)^8) &
%   arg_y := ls*(-(yaw-3.14159) + 0.16667*(yaw-3.14159)^3 - 0.0083333*(yaw-3.14159)^5 + 0.00019841*(yaw-3.14159)^7 - 2.7557e-6*(yaw-3.14159)^9) &
%   timer := 0
reset = struct('f', @Ranger_nonlinear_working_C2_St2_Tra1_ResetEq,'hasInput',1,'inputDim',3);

%% equation:
%   
A = ...
zeros([0,0]);
b = ...
zeros([0,1]);
polyOpt = struct('A', A, 'b', b);
guard = mptPolytope(polyOpt);

trans{1} = transition(guard, reset, 3,[], 1);

loc{2} = location('S2', inv, trans, dynamics);



%-------------------------State SetLinearSpeed-----------------------------

%% equation:
%   arg_x' == 0 &
%   arg_y' == 0 &
%   timer' == 1
dynA = ...
[0,0,0;0,0,0;0,0,0];
dynB = ...
[0,0,0;0,0,0;0,0,0];
dync = ...
[0;0;1];
dynamics = linearSys(dynA, dynB, dync);

%% equation:
%   timer <= 0
A = ...
[0,0,1];
b = ...
[0];
polyOpt = struct('A', A, 'b', b);
inv = mptPolytope(polyOpt);

trans = {};
%% equation:
%   timer := 0
resetA = ...
[1,0,0;0,1,0;0,0,0];
resetc = ...
[0;0;0];
reset = struct('A', resetA, 'c', resetc,'hasInput',0);

%% equation:
%   
A = ...
zeros([0,0]);
b = ...
zeros([0,1]);
polyOpt = struct('A', A, 'b', b);
guard = mptPolytope(polyOpt);

trans{1} = transition(guard, reset, 4,"set_vel", 3);

loc{3} = location('S3', inv, trans, dynamics);



%-----------------------State SetAngularSpeedArgs--------------------------

%% equation:
%   arg_x' == 0 &
%   arg_y' == 0 &
%   timer' == 1
dynA = ...
[0,0,0;0,0,0;0,0,0];
dynB = ...
[0,0,0;0,0,0;0,0,0];
dync = ...
[0;0;1];
dynamics = linearSys(dynA, dynB, dync);

%% equation:
%   timer <= 0
A = ...
[0,0,1];
b = ...
[0];
polyOpt = struct('A', A, 'b', b);
inv = mptPolytope(polyOpt);

trans = {};
%% equation:
%   arg_x := as &
%   timer := 0
resetA = ...
[0,0,0;0,1,0;0,0,0];
resetc = ...
[0;0;0];
reset = struct('A', resetA, 'c', resetc,'hasInput',0);

%% equation:
%   timer >= 0
c = [0;0;1];
d = 0;

guard = conHyperplane(c,d);

trans{1} = transition(guard, reset, 5,[]);

loc{4} = location('S4', inv, trans, dynamics);



%-------------------------State SetAngularSpeed----------------------------

%% equation:
%   arg_x' == 0 &
%   arg_y' == 0 &
%   timer' == 1
dynA = ...
[0,0,0;0,0,0;0,0,0];
dynB = ...
[0,0,0;0,0,0;0,0,0];
dync = ...
[0;0;1];
dynamics = linearSys(dynA, dynB, dync);

%% equation:
%   timer <= 0
A = ...
[0,0,1];
b = ...
[0];
polyOpt = struct('A', A, 'b', b);
inv = mptPolytope(polyOpt);

trans = {};
%% equation:
%   no reset equation given
resetA = ...
[1,0,0;0,1,0;0,0,1];
resetc = ...
[0;0;0];
reset = struct('A', resetA, 'c', resetc,'hasInput',0);

%% equation:
%   
A = ...
zeros([0,0]);
b = ...
zeros([0,1]);
polyOpt = struct('A', A, 'b', b);
guard = mptPolytope(polyOpt);

trans{1} = transition(guard, reset, 1,"set_yawVel", 3);

loc{5} = location('S5', inv, trans, dynamics);



comp{2} = hybridAutomaton(loc);

iBinds{2} = [[3,7];[1,1];[1,2]];

%---Component RoboWorld.EnvironmentMapping.Environment.RobotMovement-------

%----------------------------State Movement--------------------------------

%% equation:
%   pos_x' == vel_x & pos_y' == vel_y &
%   vel_x' == acc_x & vel_y' == acc_y &
%   acc_x' == 0 & acc_y' == 0 &
%   yaw' == yawVel & yawVel' == yawAcc & yawAcc' == 0 &
%   t' == 1
dynA = ...
[0,0,1,0,0,0,0,0,0,0;0,0,0,1,0,0,0,0,0,0;0,0,0,0,1,0,0,0,0,0;0,0,0,0,0,...
1,0,0,0,0;0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,1,0,0;0,...
0,0,0,0,0,0,0,1,0;0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0];
dynB = ...
[0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0];
dync = ...
[0;0;0;0;0;0;0;0;0;1];
dynamics = linearSys(dynA, dynB, dync);

%% equation:
%   -0.00001 <= yaw & yaw <= 6.2832
A = ...
[0,0,0,0,0,0,-1,0,0,0;0,0,0,0,0,0,1,0,0,0];
b = ...
[1.0000000000000000818030539140313E-05;6.2831999999999998962607605790254];
polyOpt = struct('A', A, 'b', b);
inv = mptPolytope(polyOpt);

trans = {};
%% equation:
%   vel_x := 0 & vel_y := 0 &
%   acc_x := 0 & acc_y := 0
resetA = ...
[1,0,0,0,0,0,0,0,0,0;0,1,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,...
0,0,0,0,0;0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,1,0,0,0;0,...
0,0,0,0,0,0,1,0,0;0,0,0,0,0,0,0,0,1,0;0,0,0,0,0,0,0,0,0,1];
resetc = ...
[0;0;0;0;0;0;0;0;0;0];
reset = struct('A', resetA, 'c', resetc,'hasInput',0);

%% equation:
%   
A = ...
zeros([0,0]);
b = ...
zeros([0,1]);
polyOpt = struct('A', A, 'b', b);
guard = mptPolytope(polyOpt);

trans{1} = transition(guard, reset, 1,"movement_blocked", 10);

%% equation:
%   vel_x := arg_x & vel_y := arg_y
resetA = ...
[1,0,0,0,0,0,0,0,0,0;0,1,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,...
0,0,0,0,0;0,0,0,0,1,0,0,0,0,0;0,0,0,0,0,1,0,0,0,0;0,0,0,0,0,0,1,0,0,0;0,...
0,0,0,0,0,0,1,0,0;0,0,0,0,0,0,0,0,1,0;0,0,0,0,0,0,0,0,0,1];
resetB = ...
[0,0;0,0;1,0;0,1;0,0;0,0;0,0;0,0;0,0;0,0];
resetc = ...
[0;0;0;0;0;0;0;0;0;0];
resetInputDim = 2;
reset = struct('A', resetA, 'B', resetB,'c', resetc,'hasInput',1,'inputDim',resetInputDim);

%% equation:
%   
A = ...
zeros([0,0]);
b = ...
zeros([0,1]);
polyOpt = struct('A', A, 'b', b);
guard = mptPolytope(polyOpt);

trans{2} = transition(guard, reset, 1,"set_vel", 10);

%% equation:
%   acc_x := arg_x & acc_y := arg_y
resetA = ...
[1,0,0,0,0,0,0,0,0,0;0,1,0,0,0,0,0,0,0,0;0,0,1,0,0,0,0,0,0,0;0,0,0,1,0,...
0,0,0,0,0;0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,1,0,0,0;0,...
0,0,0,0,0,0,1,0,0;0,0,0,0,0,0,0,0,1,0;0,0,0,0,0,0,0,0,0,1];
resetB = ...
[0,0;0,0;0,0;0,0;1,0;0,1;0,0;0,0;0,0;0,0];
resetc = ...
[0;0;0;0;0;0;0;0;0;0];
resetInputDim = 2;
reset = struct('A', resetA, 'B', resetB,'c', resetc,'hasInput',1,'inputDim',resetInputDim);

%% equation:
%   
A = ...
zeros([0,0]);
b = ...
zeros([0,1]);
polyOpt = struct('A', A, 'b', b);
guard = mptPolytope(polyOpt);

trans{3} = transition(guard, reset, 1,"set_acc", 10);

%% equation:
%   yawVel := arg_x
resetA = ...
[1,0,0,0,0,0,0,0,0,0;0,1,0,0,0,0,0,0,0,0;0,0,1,0,0,0,0,0,0,0;0,0,0,1,0,...
0,0,0,0,0;0,0,0,0,1,0,0,0,0,0;0,0,0,0,0,1,0,0,0,0;0,0,0,0,0,0,1,0,0,0;0,...
0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,1,0;0,0,0,0,0,0,0,0,0,1];
resetc = ...
[0;0;0;0;0;0;0;0;0;0];
reset = struct('A', resetA, 'c', resetc,'hasInput',0);

%% equation:
%   
A = ...
zeros([0,0]);
b = ...
zeros([0,1]);
polyOpt = struct('A', A, 'b', b);
guard = mptPolytope(polyOpt);

trans{4} = transition(guard, reset, 1,"set_yawVel", 10);

%% equation:
%   yawAcc := arg_x
resetA = ...
[1,0,0,0,0,0,0,0,0,0;0,1,0,0,0,0,0,0,0,0;0,0,1,0,0,0,0,0,0,0;0,0,0,1,0,...
0,0,0,0,0;0,0,0,0,1,0,0,0,0,0;0,0,0,0,0,1,0,0,0,0;0,0,0,0,0,0,1,0,0,0;0,...
0,0,0,0,0,0,1,0,0;0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,1];
resetc = ...
[0;0;0;0;0;0;0;0;0;0];
reset = struct('A', resetA, 'c', resetc,'hasInput',0);

%% equation:
%   
A = ...
zeros([0,0]);
b = ...
zeros([0,1]);
polyOpt = struct('A', A, 'b', b);
guard = mptPolytope(polyOpt);

trans{5} = transition(guard, reset, 1,"set_yawAcc", 10);

%% equation:
%   yaw := yaw - 6.28319
resetA = ...
[1,0,0,0,0,0,0,0,0,0;0,1,0,0,0,0,0,0,0,0;0,0,1,0,0,0,0,0,0,0;0,0,0,1,0,...
0,0,0,0,0;0,0,0,0,1,0,0,0,0,0;0,0,0,0,0,1,0,0,0,0;0,0,0,0,0,0,1,0,0,0;0,...
0,0,0,0,0,0,1,0,0;0,0,0,0,0,0,0,0,1,0;0,0,0,0,0,0,0,0,0,1];
resetc = ...
[0;0;0;0;0;0;-6.2831900000000002748379301920068;0;0;0];
reset = struct('A', resetA, 'c', resetc,'hasInput',0);

%% equation:
%   yaw >= 6.2832
c = [0;0;0;0;0;0;1;0;0;0];
d = 6.2832;

guard = conHyperplane(c,d);

trans{6} = transition(guard, reset, 1,[]);

%% equation:
%   yaw := yaw + 6.28319
resetA = ...
[1,0,0,0,0,0,0,0,0,0;0,1,0,0,0,0,0,0,0,0;0,0,1,0,0,0,0,0,0,0;0,0,0,1,0,...
0,0,0,0,0;0,0,0,0,1,0,0,0,0,0;0,0,0,0,0,1,0,0,0,0;0,0,0,0,0,0,1,0,0,0;0,...
0,0,0,0,0,0,1,0,0;0,0,0,0,0,0,0,0,1,0;0,0,0,0,0,0,0,0,0,1];
resetc = ...
[0;0;0;0;0;0;6.2831900000000002748379301920068;0;0;0];
reset = struct('A', resetA, 'c', resetc,'hasInput',0);

%% equation:
%   yaw <= -0.00001
c = [0;0;0;0;0;0;-1;0;0;0];
d = 1e-05;

guard = conHyperplane(c,d);

trans{7} = transition(guard, reset, 1,[]);

loc{1} = location('S1', inv, trans, dynamics);



comp{3} = hybridAutomaton(loc);

iBinds{3} = [[2,1];[2,2]];

%-Component RoboWorld.EnvironmentMapping.Environment.CollisionDetection----

%---------------------State OutsideCollisionRadius-------------------------

%% equation:
%   
dynA = ...
[0,0;0,0];
dynB = ...
[0;0];
dync = ...
[0;0];
dynamics = linearSys(dynA, dynB, dync);

%% equation:
%   (pos_x - obstacle_x)^2 + (pos_y - obstacle_y)^2 >= collisionRadius^2
vars = sym('x',[2,1]);
syms x1 x2;
eq = 1 - (x2 - 5)^2 - (x1 - 5)^2;
compOp = '<=';

inv = levelSet(eq,vars,compOp);

trans = {};
%% equation:
%   no reset equation given
resetA = ...
[1,0;0,1];
resetc = ...
[0;0];
reset = struct('A', resetA, 'c', resetc,'hasInput',0);

%% equation:
%   (pos_x - obstacle_x)^2 + (pos_y - obstacle_y)^2 <= collisionRadius^2
vars = sym('x',[2,1]);
syms x1 x2;
eq = (x1 - 5)^2 + (x2 - 5)^2 - 1;
compOp = '<=';

guard = levelSet(eq,vars,compOp);

trans{1} = transition(guard, reset, 1,"movement_blocked");

loc{1} = location('S1', inv, trans, dynamics);



comp{4} = hybridAutomaton(loc);

iBinds{4} = [0 1];
%Component RoboWorld.EnvironmentMapping.Environment.Obstacle_InputEventMapping

%------------------------State InDetectionRadius---------------------------

%% equation:
%   
dynA = ...
[0,0;0,0];
dynB = ...
[0;0];
dync = ...
[0;0];
dynamics = linearSys(dynA, dynB, dync);

%% equation:
%   (pos_x - obstacle_x)^2 + (pos_y - obstacle_y)^2 <= (detectionRadius + resetRadius)^2
vars = sym('x',[2,1]);
syms x1 x2;
eq = (x1 - 5)^2 + (x2 - 5)^2 - 1266645839461327/140737488355328;
compOp = '<=';

inv = levelSet(eq,vars,compOp);

trans = {};
%% equation:
%   no reset equation given
resetA = ...
[1,0;0,1];
resetc = ...
[0;0];
reset = struct('A', resetA, 'c', resetc,'hasInput',0);

%% equation:
%   (pos_x - obstacle_x)^2 + (pos_y - obstacle_y)^2 >= (detectionRadius + resetRadius)^2
vars = sym('x',[2,1]);
syms x1 x2;
eq = 1266645839461327/140737488355328 - (x2 - 5)^2 - (x1 - 5)^2;
compOp = '<=';

guard = levelSet(eq,vars,compOp);

trans{1} = transition(guard, reset, 2,[]);

loc{1} = location('S1', inv, trans, dynamics);



%---------------------State OutsideDetectionRadius-------------------------

%% equation:
%   
dynA = ...
[0,0;0,0];
dynB = ...
[0;0];
dync = ...
[0;0];
dynamics = linearSys(dynA, dynB, dync);

%% equation:
%   (pos_x - obstacle_x)^2 + (pos_y - obstacle_y)^2 >= detectionRadius^2
vars = sym('x',[2,1]);
syms x1 x2;
eq = 9 - (x2 - 5)^2 - (x1 - 5)^2;
compOp = '<=';

inv = levelSet(eq,vars,compOp);

trans = {};
%% equation:
%   no reset equation given
resetA = ...
[1,0;0,1];
resetc = ...
[0;0];
reset = struct('A', resetA, 'c', resetc,'hasInput',0);

%% equation:
%   (pos_x - obstacle_x)^2 + (pos_y - obstacle_y)^2 <= detectionRadius^2
vars = sym('x',[2,1]);
syms x1 x2;
eq = (x1 - 5)^2 + (x2 - 5)^2 - 9;
compOp = '<=';

guard = levelSet(eq,vars,compOp);

trans{1} = transition(guard, reset, 1,"obstacle");

loc{2} = location('S2', inv, trans, dynamics);



comp{5} = hybridAutomaton(loc);

iBinds{5} = [0 1];
HA = parallelHybridAutomaton(comp,iBinds);


end