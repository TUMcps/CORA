function display(sys)
% display - Displays a linProbSys object on the command window
%
% Syntax:
%    display(sys)
%
% Inputs:
%    sys - linProbSys object
%
% Outputs:
%    -
%
% Example:
%    A = [-1 -4; 4 -1];
%    B = eye(2);
%    C = 0.7*eye(2);
%    sys = linProbSys(A,B,C)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       06-October-2007
% Last update:   19-June-2022
%                23-November-2022 (TL, dispInput)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% disp input if necessary
dispInput(inputname(1))

%display parent object
display@contDynamics(sys);

%display type
disp('Type: Linear continuous-time probabilistic system');

% display state matrix
disp("System matrix:");
displayMatrixVector(sys.A,"A");

% display input matrix
disp("Input matrix:");
displayMatrixVector(sys.B,"B");

% display noise matrix offset
disp("Noise matrix:");
displayMatrixVector(sys.C,"C");

% ------------------------------ END OF CODE ------------------------------
