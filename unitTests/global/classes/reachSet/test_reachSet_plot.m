function res = test_reachSet_plot
% test_reachSet_plot - unit test function for plot
%
% Syntax:
%    res = test_reachSet_plot()
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

% Authors:       Mark Wetzlinger
% Written:       01-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init result
res = true;

figure;

% empty reachSet
R = reachSet();
try
    plot(R);
    res = false;
end

% discrete-time system
A = [0.9810    0.0143    0.0262   -0.0140;
    0.0079    1.0133    0.0267    0.0416;
    0.0029   -0.0054    0.9680    0.0188;
    0.0503   -0.0242    0.0411    0.9877];
B = [-0.0246    0.0003    0.0372;
    0.0251    0.0121    0.0006;
   -0.0001   -0.0003    0.0125;
   -0.0009   -0.0250    0.0136];
dt = 0.05;
sys = linearSysDT(A,B,dt);
params.R0 = zonotope([10;-5;8;-12],0.5*eye(4));
params.tFinal = 5;
options.zonotopeOrder = 50;

% compute (only time-point!) reachable set
R = reach(sys,params,options);

try
    % plot
    plot(R);
    plot(R,[2,3]);

    % with linespec
    plot(R,[1,2],'r');

    % with name-value pairs
    plot(R,[1,2],'Order',10);

catch
    res = false;
end

% linear continuous-time system
A = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
sys = linearSys(A,1);
params.R0 = zonotope(ones(5,1),0.1*eye(5));
params.U = zonotope([1; 0; 0; 0.5; -0.5],0.5*diag([0.2, 0.5, 0.2, 0.5, 0.5]));
params.tFinal = 2;
options.linAlg = 'adaptive';
options.error = 0.1;

% compute reachable set
R = reach(sys,params,options);

try
    % plot
    plot(R);
    plot(R,[1,2]);

    % with linespec
    plot(R,[1,2],'r');

    % with name-value pairs
    plot(R,[1,2],'Order',10);

    % explicitly, time-point solution
    plot(R,[1,2],'Set','tp');

catch
    res = false;
end

% close figure
close

% ------------------------------ END OF CODE ------------------------------
