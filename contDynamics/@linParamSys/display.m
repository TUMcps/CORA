function display(sys)
% display - Displays a linParamSys object on the command window
%
% Syntax:
%    display(sys)
%
% Inputs:
%    sys - linParamSys object
%
% Outputs:
%    ---
%
% Example:
%    Ac = [-2 0; 1.5 -3];
%    Aw = [0 0; 0.5 0];
%    A = intervalMatrix(Ac,Aw);
%    B = [1; 1];
%    sys = linParamSys(A,B,'varParam')
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       05-August-2010
% Last update:   19-June-2022
%                23-November-2022 (TL, dispInput)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% disp input if necessary
dispInput(inputname(1))

%display parent object
display@contDynamics(sys);

%display type
disp('Type: Linear continuous-time parametric system');

% display parameter type
if sys.constParam
    disp('  parameter: constant'); 
else
    disp('  parameter: time-varying');
end

% display state matrix
disp("System matrix:");
displayMatrixVector(sys.A,"A");

% display input matrix
disp("Input matrix:");
displayMatrixVector(sys.B,"B");

% ------------------------------ END OF CODE ------------------------------
