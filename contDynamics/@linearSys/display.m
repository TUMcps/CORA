function display(linsys)
% display - Displays a linearSys object on the command window
%
% Syntax:
%    display(linsys)
%
% Inputs:
%    linsys - linearSys object
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
%    linsys = linearSys(A,B,c,C,D,k)
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
display@contDynamics(linsys);

%display type
disp("Type: Linear continuous-time time-invariant system");

% state equation
disp("State-space equation: x' = Ax + Bu + c + Ew");

% display state matrix
disp("System matrix:");
displayMatrixVector(linsys.A,"A");

% display input matrix
disp("Input matrix:");
displayMatrixVector(linsys.B,"B");

% display constant offset
disp("Constant offset:");
displayMatrixVector(linsys.c,"c");

% display constant offset
disp("Disturbance matrix:");
displayMatrixVector(linsys.E,"E");

% check if there is an output equation
isOutput = ~isscalar(linsys.C) || linsys.C ~= 1 || any(any(linsys.D)) ...
    || any(linsys.k) || any(any(linsys.F));

% output equation
if ~isOutput
    disp("Output equation: y = x");
else
    disp("Output equation: y = Cx + Du + k + Fv");
    
    % display output matrix
    disp("Output matrix:");
    displayMatrixVector(linsys.C,"C");
    
    % display feedthrough matrix
    disp("Feedthrough matrix:");
    displayMatrixVector(linsys.D,"D");
    
    % display constant offset
    disp("Constant offset:");
    displayMatrixVector(linsys.k,"k");

    % display noise matrix
    disp("Noise matrix:");
    displayMatrixVector(linsys.F,"F");
end

% ------------------------------ END OF CODE ------------------------------
