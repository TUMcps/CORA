function display(sys)
% display - displays a contDynamics object on the command window
%
% Syntax:
%    display(sys)
%
% Inputs:
%    sys - contDynamics object
%
% Outputs:
%    ---
%
% Example: 
%    sys = contDynamics('test',3,1,2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       02-May-2007
% Last update:   17-October-2007
%                19-June-2022
%                23-November-2022 (TL, dispInput)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% disp input if necessary
dispInput(inputname(1))

% display name and id
disp("Continuous dynamics: '" + sys.name + "'");
% display dynamic properties
disp("  number of dimensions: " + sys.nrOfDims);
disp("  number of inputs: " + sys.nrOfInputs);
disp("  number of outputs: " + sys.nrOfOutputs);
disp("  number of disturbances: " + sys.nrOfDisturbances);
disp("  number of noises: " + sys.nrOfNoises);

fprintf(newline);

% ------------------------------ END OF CODE ------------------------------
