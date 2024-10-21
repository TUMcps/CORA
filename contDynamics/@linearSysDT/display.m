function display(linsysDT)
% display - Displays a linearSysDT object on the command window
%
% Syntax:
%    display(linsysDT)
%
% Inputs:
%    linsysDT - linearSysDT object
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
%    linsys = linearSys(A,B,c,C,D,k);
% 
%    % discretize continuous-time system
%    dt = 0.05;
%    linsysDT = linearSysDT(linsys,dt);
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
display@contDynamics(linsysDT);

%display type
disp("Type: Linear discrete-time system");

% display sampling time
disp("Sampling time: " + linsysDT.dt);

% state equation
disp("State-space equation: x_k+1 = Ax_k + Bu_k + c + Ew_k");

% display state matrix
disp("System matrix:");
displayMatrixVector(linsysDT.A,"A");

% display input matrix
disp("Input matrix:");
displayMatrixVector(linsysDT.B,"B");

% display constant offset
disp("Constant offset:");
displayMatrixVector(linsysDT.c,"c");

% display constant offset
disp("Disturbance matrix:");
displayMatrixVector(linsysDT.E,"E");

% check if there is an output equation
isOutput = ~isscalar(linsysDT.C) || linsysDT.C ~= 1 || any(any(linsysDT.D)) ...
    || any(linsysDT.k) || any(any(linsysDT.F));

% output equation
if ~isOutput
    disp("Output equation: y_k = x_k");
else
    disp("Output equation: y_k = Cx_k + Du_k + k + Fv_k");
    
    % display output matrix
    disp("Output matrix:");
    displayMatrixVector(linsysDT.C,"C");
    
    % display feedthrough matrix
    disp("Feedthrough matrix:");
    displayMatrixVector(linsysDT.D,"D");
    
    % display constant offset
    disp("Constant offset:");
    displayMatrixVector(linsysDT.k,"k");

    % display noise matrix
    disp("Noise matrix:");
    displayMatrixVector(linsysDT.F,"F");
end

% ------------------------------ END OF CODE ------------------------------
