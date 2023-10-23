function display(sys)
% display - Displays a linearSysDT object on the command window
%
% Syntax:
%    display(sys)
%
% Inputs:
%    sys - linearSysDT object
%
% Outputs:
%    ---
%
% Example:
%    A = [-0.3780    0.2839    0.5403   -0.2962
%          0.1362    0.2742    0.5195    0.8266
%          0.0502   -0.1051   -0.6572    0.3874
%          1.0227   -0.4877    0.8342   -0.2372];
%    B = 0.25 * [-2 0 3; 2 1 0; 0 0 1; 0 -2 1];
%    c = 0.05 * [-4; 2; 3; 1];
%    C = [1 1 0 0; 0 -0.5 0.5 0];
%    D = [0 0 1; 0 0 0];
%    k = [0; 0.02];
%    sys_CT = linearSys(A,B,c,C,D,k);
% 
%    % discretize continuous-time system
%    dt = 0.05;
%    sys = linearSysDT(sys_CT,dt);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       16-May-2007
% Last update:   19-June-2022
%                23-November-2022 (TL, dispInput)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% disp input if necessary
dispInput(inputname(1))

%display parent object
display@contDynamics(sys);

%display type
disp("Type: Linear discrete-time system");

% display sampling time
disp("Sampling time: " + sys.dt);

% state equation
disp("x' = Ax + Bu + c");

% display state matrix
disp("System matrix:");
displayMatrixVector(sys.A,"A");

% display input matrix
disp("Input matrix:");
displayMatrixVector(sys.B,"B");

% display constant offset
disp("Constant offset:");
displayMatrixVector(sys.c,"c");

% check if there is an output equation
isOutput = ~isempty(sys.C) || ~isempty(sys.D) || ~isempty(sys.k);

% output equation
if isOutput
    disp("y = Cx + Du + k");
    
    % display output matrix
    disp("Output matrix:");
    displayMatrixVector(sys.C,"C");
    
    % display feedthrough matrix
    disp("Feedthrough matrix:");
    displayMatrixVector(sys.D,"D");
    
    % display constant offset
    disp("Constant offset:");
    displayMatrixVector(sys.k,"k");
end

% ------------------------------ END OF CODE ------------------------------
